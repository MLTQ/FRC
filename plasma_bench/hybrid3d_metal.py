"""Optional fused Metal kernels, with explicit single-precision boundaries."""
import time
import numpy as np


class MetalKernels:
    """Transport fluxes and CIC moments; host integration remains float64."""
    def __init__(self, model):
        import mlx.core as mx
        if not mx.metal.is_available():
            raise RuntimeError('Metal backend requested but no GPU is available')
        self.mx, self.model = mx, model
        self.shape = tuple(model.box.shape)
        self.size = int(np.prod(self.shape))
        self.calls = self.wall_s = 0
        nx, ny, nz = self.shape
        header = f'const int nx={nx}, ny={ny}, nz={nz}, N={self.size};'
        transport = header + r'''
            uint q = thread_position_in_grid.x;
            if (q >= N) return;
            int i = q/(ny*nz), j = (q/nz)%ny, k = q%nz;
            int plus[3] = {int(q)+(i+1<nx ? ny*nz : -(nx-1)*ny*nz),
                           int(q)+(j+1<ny ? nz : -(ny-1)*nz),
                           int(q)+(k+1<nz ? 1 : -(nz-1))};
            int minus[3] = {int(q)+(i>0 ? -ny*nz : (nx-1)*ny*nz),
                            int(q)+(j>0 ? -nz : (ny-1)*nz),
                            int(q)+(k>0 ? -1 : nz-1)};
            float boundary = 0.0f, outgoing = 0.0f;
            for (int a=0; a<3; ++a) {
                float speed = 0.5f*(velocity[a*N+q]+velocity[a*N+plus[a]]);
                float upstream = speed>=0.0f ? u[q] : u[plus[a]];
                out[a*N+q] = links[a*N+q]*speed*upstream;
                float previous = 0.5f*(velocity[a*N+minus[a]]+velocity[a*N+q]);
                outgoing += (links[a*N+q]*max(speed,0.0f)
                    +links[a*N+minus[a]]*max(-previous,0.0f))*scales[3+a];
                float va = velocity[a*N+q];
                boundary += bp[a*N+q]*max(va,0.0f)+bm[a*N+q]*max(-va,0.0f);
            }
            out[3*N+q] = boundary;
            out[4*N+q] = outgoing;
        '''
        self.transport_kernel = mx.fast.metal_kernel(name='frc_energy_transport',
            input_names=['u', 'velocity', 'links', 'bp', 'bm', 'scales'], output_names=['out'],
            source=transport, compile_options={'math_mode': 'safe'})
        deposit = header + r'''
            uint p = thread_position_in_grid.x;
            if (p >= positions_shape[0]) return;
            float sx=(positions[3*p]-scales[0])*scales[3];
            float sy=(positions[3*p+1]-scales[1])*scales[4];
            float sz=(positions[3*p+2]-scales[2])*scales[5];
            int ix=int(floor(sx)), iy=int(floor(sy)), iz=int(floor(sz));
            float ax=sx-floor(sx), ay=sy-floor(sy), az=sz-floor(sz);
            for (int c=0;c<8;++c) {
                int i=((ix+((c&1)!=0))%nx+nx)%nx;
                int j=((iy+((c&2)!=0))%ny+ny)%ny;
                int k=((iz+((c&4)!=0))%nz+nz)%nz;
                int q=(i*ny+j)*nz+k;
                float w=weights[p]*scales[6]*(c&1 ? ax : 1.0f-ax)
                    *(c&2 ? ay : 1.0f-ay)*(c&4 ? az : 1.0f-az);
                // Combine equal-node contributions in groups of eight lanes.
                // This is only an arithmetic optimization: unequal keys use
                // the general atomic path, regardless of particle ordering.
                uint lane=thread_index_in_simdgroup;
                bool combine=false;
                if ((p/8)*8+7<positions_shape[0]) {
                    uint leader=simd_shuffle(uint(q),lane&~7u);
                    uint same=uint(uint(q)==leader);
                    same+=simd_shuffle_xor(same,1);
                    same+=simd_shuffle_xor(same,2);
                    same+=simd_shuffle_xor(same,4);
                    combine=(same==8);
                }
                for (int a=0;a<4;++a) {
                    float value=a==0 ? w : w*velocities[3*p+a-1]*scales[7];
                    if (combine) {
                        value+=simd_shuffle_xor(value,1);
                        value+=simd_shuffle_xor(value,2);
                        value+=simd_shuffle_xor(value,4);
                        if ((lane&7u)==0)
                            atomic_fetch_add_explicit(&moments[a*N+q],value,memory_order_relaxed);
                    } else {
                        atomic_fetch_add_explicit(&moments[a*N+q],value,memory_order_relaxed);
                    }
                }
            }
        '''
        self.deposit_kernel = mx.fast.metal_kernel(name='frc_cic_moments',
            input_names=['positions', 'velocities', 'weights', 'scales'], output_names=['moments'],
            source=deposit, atomic_outputs=True, compile_options={'math_mode': 'safe'})
        self.density_scale = max(float(model.n_ref), 1.)
        self.speed_scale = 1e5
        self.scales = self.array([*model.lo, *list(1/model.d),
            1/(self.density_scale*model.box.cell_volume), 1/self.speed_scale])

    def array(self, value):
        return self.mx.array(np.ascontiguousarray(value, dtype=np.float32))

    def prepare_geometry(self, electrons):
        self.geometry = [self.array(v) for v in
                         (electrons.link_area, electrons.boundary_plus, electrons.boundary_minus)]

    def transport(self, u, velocity):
        """Return single-precision face fluxes and wall rate, synchronously."""
        start = time.perf_counter()
        result = self.transport_kernel(inputs=[self.array(u), self.array(velocity), *self.geometry, self.scales],
            output_shapes=[(5,)+self.shape], output_dtypes=[self.mx.float32],
            grid=(self.size, 1, 1), threadgroup=(min(256, self.size), 1, 1))[0]
        self.mx.eval(result)
        value = np.asarray(result).astype(np.float64)
        self.calls += 1; self.wall_s += time.perf_counter()-start
        return value[:3], value[3], value[4]*(1+1e-6)

    def prepare_moments(self, x, v, w):
        return [self.array(x), self.array(v), self.array(w), self.scales]

    def moments_device(self, inputs):
        """Device-only work, separately benchmarked from uploads/collection."""
        count = int(inputs[0].shape[0])
        return self.deposit_kernel(inputs=inputs,
            output_shapes=[(4,)+self.shape], output_dtypes=[self.mx.float32],
            grid=(count, 1, 1), threadgroup=(min(256, count), 1, 1), init_value=0)[0]

    def moments(self, x, v, w):
        """Scaled float32 atomics; position integration and result are float64."""
        if len(x) == 0:
            return np.zeros((4,)+self.shape)
        start = time.perf_counter()
        result = self.moments_device(self.prepare_moments(x, v, w))
        self.mx.eval(result)
        value = np.asarray(result).astype(np.float64)
        value *= self.density_scale; value[1:] *= self.speed_scale
        self.calls += 1; self.wall_s += time.perf_counter()-start
        return value
