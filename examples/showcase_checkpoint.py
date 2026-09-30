"""Render the latest frozen physical fields without advancing the simulation."""
import argparse,base64,hashlib,json
from pathlib import Path
from functools import partial
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .regular_coupling_reference import initialize
from .compatible_field_norms import direct_basis


def build():
    root=Path(__file__).resolve().parents[1];out=root/'outputs/frc_showcase_checkpoint';out.mkdir(exist_ok=True)
    source=root/'outputs/regular_coupling_verification/case_07.npz';record=json.loads(source.with_suffix('.json').read_text());digest=hashlib.sha256(source.read_bytes()).hexdigest();assert digest==record['arrays_sha256']
    a=np.load(source);ref,d,u,m,field=initialize(record['spec']);mesh=m.mesh
    nr,nz=128,256;r=(np.arange(nr)+.5)*ref.radius/nr;z=(np.arange(nz)+.5)*2*ref.half_length/nz-ref.half_length
    nx,dx,ex=direct_basis(.5*mesh.rf**2,.5*r*r,5);zz,dz,ez=direct_basis(mesh.zf,z,5,True);nodes=(mesh.nr+1)*mesh.nz;nt=(mesh.nr-1)*mesh.nz
    j=a['current'];jr=r[:,None]*(ex@j[:mesh.nr*mesh.nz].reshape(mesh.nr,mesh.nz)@zz.T)
    theta=np.pad(j[mesh.nr*mesh.nz:mesh.nr*mesh.nz+nt].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
    jt=(nx@(mesh.rf[:,None]*theta)@zz.T)/r[:,None];jz=nx@j[mesh.nr*mesh.nz+nt:].reshape(mesh.nr+1,mesh.nz)@ez.T
    magnetic=a['magnetic'];p=magnetic[:nodes].reshape(mesh.nr+1,mesh.nz);gg=magnetic[nodes:].reshape(mesh.nr,mesh.nz)
    b=np.stack((-(nx@p@dz.T)/r[:,None],r[:,None]*(ex@gg@ez.T),dx@p@zz.T+ref.axial_bias),axis=-1)
    jm=np.sqrt(jr*jr+jt*jt+jz*jz);bm=np.linalg.norm(b,axis=-1);scale=1/np.hypot(ref.radius,ref.half_length)
    interp=RegularGridInterpolator((r,z),b,bounds_error=False,fill_value=None)
    paths=[];step=2*ref.half_length/230
    for rad in ref.radius*np.array([.12,.3,.5,.7,.9]):
        for angle in np.arange(6)*2*np.pi/6:
            value=np.array([rad,angle,z[0]]);points=[];colors=[]
            def velocity(q):
                bb=interp([[q[0],q[2]]])[0];norm=np.linalg.norm(bb)
                return np.array([bb[0],bb[1]/q[0],bb[2]])/norm
            for _ in range(260):
                if not(r[0]<=value[0]<=r[-1] and z[0]-1e-14<=value[2]<=z[-1]):break
                rr,phi,ax=value;points.append((scale*np.array([ax,rr*np.cos(phi),rr*np.sin(phi)])).tolist());colors.append(float(np.linalg.norm(interp([[rr,ax]])[0])*1000))
                k1=velocity(value);k2=velocity(value+.5*step*k1);k3=velocity(value+.5*step*k2);k4=velocity(value+step*k3);value+=step*(k1+2*k2+2*k3+k4)/6
            paths.append(dict(points=points,b=colors))
    texture=np.stack((jm,bm,nx@p@zz.T,b[...,0]),axis=-1).transpose(1,0,2).astype('<f4')
    data=dict(source=str(source.relative_to(root)),source_sha256=digest,calculation_source_sha256=record['source_sha256'],texture_size=[nr,nz],extent=[ref.radius*scale,ref.half_length*scale],jmax=float(jm.max()),bmax_mT=float(bm.max()*1000),bmin_mT=float(bm.min()*1000),frames=[dict(time_ms=0,bias_mT=ref.axial_bias*1000,texture=base64.b64encode(texture.tobytes()).decode(),rings=[],lines=[],lines3d=paths)],interpretation='Frozen instantaneous verification fixture; no new physical evolution or formation claim. Animation is camera orbit.',physical_radius_m=ref.radius,physical_length_m=2*ref.half_length,uniform_display_scale_per_m=scale)
    (out/'data.json').write_text(json.dumps(data,separators=(',',':')))
    assets=root/'visualizations/frc_showcase';renderer=(assets/'renderer.js').read_text().replace('this.setFrame(2);','this.setFrame(Math.min(2,data.frames.length-1));')
    renderer=renderer.replace('vec2 t=vec2(length(p.yz)/extent.x,','if(length(p.yz)>extent.x||abs(p.x)>extent.y)continue;vec2 t=vec2(length(p.yz)/extent.x,')
    renderer=renderer.replace('else {for(const path of f.lines)',"else if(f.lines3d){for(const path of f.lines3d)tube(vertices,path.points,path.b.map(b=>(b-data.bmin_mT)/(data.bmax_mT-data.bmin_mT)));}\nelse {for(const path of f.lines)")
    controller=(assets/'controller.js').read_text().replace('distance:5.55','distance:3.65').replace('frame:2,orbit:false','frame:0,orbit:true')
    replacements={'Inside the field reversal.':'Inside the verified fields.','FIELD STUDY 01':'RESEARCH CHECKPOINT','Current and magnetic structure from an initialized FRC simulation':'Latest verified current and magnetic fields · frozen instantaneous fixture','Toroidal current volume · color = |J|':'All current components · volume color = |J|','Poloidal field lines · color = |B|':'Full magnetic vector · color = |B|','INITIAL RADIUS  1 m':'RADIUS  35 mm  /  LENGTH  240 mm','Initialized, axisymmetric verification case — not demonstrated formation or 3D stability.':'Verified field snapshot. Camera animation only; no new formation or stability claim.','Fixed color scales. Illustrative opacity and field-line spacing. Heated-edge accuracy remains unresolved.':'402 regression tests + 16 audited coupling cases. Current transport candidate: 4 focused tests; full audit pending.','ARCHIVED CASE  source_g1_32  /  32 × 64 MATERIAL CELLS  /  EXACTLY ZERO EXTERIOR PLASMA CURRENT':'REGULAR COUPLING CASE 07 / 32 × 32 FIELD CELLS / FULL TOROIDAL + POLOIDAL MAGNETIC COMPONENTS','outputs/frc_showcase/poster.png':'outputs/frc_showcase_checkpoint/poster.png'}
    for old,new in replacements.items():controller=controller.replace(old,new)
    controller=controller.replace('(f*max).toFixed(2)',"(f*(id==='b-ticks'?DATA.bmax_mT-DATA.bmin_mT:max)+(id==='b-ticks'?DATA.bmin_mT:0)).toFixed(2)")
    controller=controller.replace("(k/4*(i?DATA.bmax_mT:DATA.jmax/1000)).toFixed(2)","(k/4*(i?DATA.bmax_mT-DATA.bmin_mT:DATA.jmax/1000)+(i?DATA.bmin_mT:0)).toFixed(2)")
    html=(assets/'index.html').read_text().replace('FRC · Inside the field reversal','FRC · Verified field checkpoint')
    for old,new in replacements.items():html=html.replace(old,new)
    html=html.replace('INITIALIZED FRC TEST','VERIFIED FIELD SNAPSHOT').replace('32 × 64 material cells','32 × 32 field cells').replace('Plasma current and magnetic structure, reconstructed from a saved simulation.','The latest verified fields, with full toroidal and poloidal magnetic structure.')
    html=html.replace('<button data-frame="0">0 ms</button><button data-frame="1">20 ms</button><button class="selected" data-frame="2">40 ms</button>','<button class="selected" data-frame="0">Snapshot</button>').replace('id="orbit" aria-pressed="false"','id="orbit" aria-pressed="true"')
    html=html.replace('Toroidal current · volume rendering of |J|','All three current components · volume rendering of |J|').replace('Poloidal field lines · colored by |B|','Full magnetic streamlines · colored by |B|').replace('CLOSED LOOPS + SURROUNDING VACUUM FIELD','FULL TOROIDAL + POLOIDAL MAGNETIC FIELD')
    html=html.replace('A plasma initialized with reversed magnetic flux, then evolved with resistive heating and a changing bias field. <strong>This is a verification experiment, not demonstrated FRC formation.</strong>','A frozen fixture from the latest independently audited fluid–magnetic coupling. <strong>402 regression tests and 16 audited cases.</strong> New current transport passes 4 focused tests; its full audit is saved for next week.')
    html=html.replace('Initial radius <strong>1 m</strong>','Radius <strong>35 mm</strong> · length <strong>240 mm</strong>').replace('40 ms','0 ms').replace('magnetic vacuum outside','filled periodic cylinder').replace('Mathematical test parameters; heated-edge accuracy remains unresolved.','Prescribed verification fields; this is not an FRC formation trajectory.').replace('Both panels share a camera. Colors retain the same scale at every saved time. The 3D shape is an axisymmetric revolution; it does not demonstrate 3D stability.','Auto orbit animates the camera around one saved state. Streamlines include the real toroidal field. Geometry is scaled uniformly; opacity and seed spacing are illustrative.')
    html=html.replace('FRC RESEARCH / ARCHIVED CASE source_g1_32 / NO EXTERIOR PLASMA','FRC RESEARCH / COUPLING CASE 07 / SOURCE 8c9339a6 / PHYSICS PAUSED FOR NEXT WEEK')
    html=html.replace('<button id="export">','<button id="video">Save orbit animation</button><button id="export">')
    video='''
document.getElementById('video').onclick=async()=>{const old=state.orbit;state.orbit=false;status.textContent='Recording 6-second camera orbit…';
const c=document.createElement('canvas');c.width=1600;c.height=900;const ctx=c.getContext('2d');const stream=c.captureStream(24);const chunks=[];const rec=new MediaRecorder(stream,{mimeType:'video/webm;codecs=vp9',videoBitsPerSecond:5000000});rec.ondataavailable=e=>chunks.push(e.data);const ended=new Promise(r=>rec.onstop=r);rec.start();
for(let i=0;i<144;i++){state.yaw=.58+i*2*Math.PI/144;ctx.fillStyle='#080e17';ctx.fillRect(0,0,1600,900);ctx.fillStyle='#e4eeee';ctx.font='30px sans-serif';ctx.fillText('VIRTUAL PLASMA BENCH · VERIFIED FIELD SNAPSHOT',45,55);ctx.font='18px sans-serif';ctx.fillStyle='#a4bbc8';ctx.fillText('CURRENT DENSITY  |J|',45,104);ctx.fillText('MAGNETIC FLUX DENSITY  |B|',830,104);
for(let k=0;k<2;k++){views[k].draw(state,[760,640]);ctx.drawImage(views[k].canvas,30+800*k,120);}ctx.fillStyle='#a4bbc8';ctx.font='19px sans-serif';ctx.fillText(`|J|: 0–${(DATA.jmax/1000).toFixed(2)} kA/m²`,45,790);ctx.fillText(`|B|: ${DATA.bmin_mT.toFixed(3)}–${DATA.bmax_mT.toFixed(3)} mT`,830,790);ctx.fillText('Radius 35 mm · length 240 mm · camera orbit only · no new plasma evolution',45,840);ctx.font='16px sans-serif';ctx.fillText('402 regression tests · 16 audited coupling cases · full moving-species transport remains open',45,876);await new Promise(r=>setTimeout(r,1000/24));}
rec.stop();await ended;stream.getTracks().forEach(t=>t.stop());const blob=new Blob(chunks,{type:'video/webm'});await fetch('/save-orbit',{method:'POST',body:blob});status.textContent='Orbit animation saved';state.orbit=old;dirty=true;};
'''
    html=html.replace('/* DATA */','const DATA='+json.dumps(data,separators=(',',':'))+';').replace('/* RENDERER */',renderer).replace('/* CONTROLLER */',controller+video)
    (out/'index.html').write_text(html);(out/'manifest.json').write_text(json.dumps({k:v for k,v in data.items() if k!='frames'},indent=2)+'\n');return out


class Handler(SimpleHTTPRequestHandler):
    def do_POST(self):
        names={'/save-poster':('poster.png',b'\x89PNG\r\n\x1a\n'),'/save-orbit':('orbit.webm',b'\x1aE\xdf\xa3')};item=names.get(self.path);n=int(self.headers.get('Content-Length','0'))
        if item is None or not 0<n<30_000_000:self.send_error(400);return
        data=self.rfile.read(n)
        if not data.startswith(item[1]):self.send_error(400);return
        (Path(self.directory)/item[0]).write_bytes(data);self.send_response(200);self.end_headers();self.wfile.write(b'Saved')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--serve',action='store_true');p.add_argument('--port',type=int,default=8779);args=p.parse_args();out=build();print(out,flush=True)
    if args.serve:ThreadingHTTPServer(('127.0.0.1',args.port),partial(Handler,directory=str(out))).serve_forever()
