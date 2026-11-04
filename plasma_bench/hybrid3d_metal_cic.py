"""Optional CIC prototype: CPU64 cell coordinates, Metal32 local weights."""
import numpy as np
from numba import njit, prange

from .hybrid3d_metal import MetalKernels


@njit(parallel=True, cache=True)
def split_coordinates(x, lo, spacing):
    cells = np.empty(x.shape, dtype=np.int32)
    fractions = np.empty(x.shape, dtype=np.float32)
    for p in prange(len(x)):
        for a in range(3):
            coordinate = (x[p, a]-lo[a])/spacing[a]
            cell = np.floor(coordinate)
            cells[p, a] = int(cell)
            fractions[p, a] = coordinate-cell
    return cells, fractions


class SplitCoordinateMetal(MetalKernels):
    """Keep absolute coordinate cancellation out of GPU32 interpolation."""
    def __init__(self, model):
        super().__init__(model)
        nx, ny, nz = self.shape
        source = f'const int nx={nx}, ny={ny}, nz={nz}, N={self.size};' + r'''
            uint p=thread_position_in_grid.x;
            if (p>=cells_shape[0]) return;
            int ix=cells[3*p], iy=cells[3*p+1], iz=cells[3*p+2];
            float ax=fractions[3*p], ay=fractions[3*p+1], az=fractions[3*p+2];
            for (int c=0;c<8;++c) {
                int i=((ix+((c&1)!=0))%nx+nx)%nx;
                int j=((iy+((c&2)!=0))%ny+ny)%ny;
                int k=((iz+((c&4)!=0))%nz+nz)%nz;
                int q=(i*ny+j)*nz+k;
                float w=weights[p]*scales[6]*(c&1 ? ax : 1.0f-ax)
                    *(c&2 ? ay : 1.0f-ay)*(c&4 ? az : 1.0f-az);
                uint lane=thread_index_in_simdgroup;
                bool combine=false;
                if ((p/8)*8+7<cells_shape[0]) {
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
        self.split_kernel = self.mx.fast.metal_kernel(name='frc_cic_split_coordinates',
            input_names=['cells','fractions','velocities','weights','scales'],
            output_names=['moments'], source=source, atomic_outputs=True,
            compile_options={'math_mode':'safe'})

    def prepare_moments(self, x, v, w):
        cells, fractions = split_coordinates(x, self.model.lo, self.model.d)
        return [self.mx.array(cells), self.mx.array(fractions),
                self.array(v), self.array(w), self.scales]

    def moments_device(self, inputs):
        count = int(inputs[0].shape[0])
        return self.split_kernel(inputs=inputs, output_shapes=[(4,)+self.shape],
            output_dtypes=[self.mx.float32], grid=(count,1,1),
            threadgroup=(min(256,count),1,1), init_value=0)[0]
