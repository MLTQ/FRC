const state={yaw:.58,pitch:.30,distance:5.55,density:.85,cutaway:true,frame:2,orbit:false};
let views=[],dirty=true,last=0;
const status=document.getElementById('status');
function redraw(){for(const v of views)v.draw(state);dirty=false;}
function metrics(){const f=DATA.frames[state.frame];document.getElementById('time-readout').textContent=f.time_ms+' ms';document.getElementById('bias-readout').textContent=f.bias_mT.toFixed(3)+' mT';}
try{
views=[new FieldView(document.getElementById('current'),DATA,0),new FieldView(document.getElementById('magnetic'),DATA,1)];
for(const [id,max] of [['j-ticks',DATA.jmax/1000],['b-ticks',DATA.bmax_mT]])document.getElementById(id).innerHTML=[0,.25,.5,.75,1].map(f=>`<span>${(f*max).toFixed(2)}</span>`).join('');metrics();
for(const canvas of document.querySelectorAll('canvas')){let drag=null;
canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);canvas.style.cursor='grabbing';};
canvas.onpointermove=e=>{if(!drag)return;state.yaw-=(e.clientX-drag[0])*.006;state.pitch=Math.max(-1.25,Math.min(1.25,state.pitch+(e.clientY-drag[1])*.006));drag=[e.clientX,e.clientY];dirty=true;};
canvas.onpointerup=()=>{drag=null;canvas.style.cursor='grab';};canvas.onwheel=e=>{e.preventDefault();state.distance=Math.max(2.5,Math.min(7,state.distance*Math.exp(e.deltaY*.001)));dirty=true;};}
document.querySelectorAll('[data-frame]').forEach(b=>b.onclick=()=>{state.frame=+b.dataset.frame;for(const v of views)v.setFrame(state.frame);document.querySelectorAll('[data-frame]').forEach(x=>x.classList.toggle('selected',x===b));metrics();dirty=true;});
document.getElementById('orbit').onclick=e=>{state.orbit=!state.orbit;e.currentTarget.setAttribute('aria-pressed',state.orbit);};
document.getElementById('cutaway').onclick=e=>{state.cutaway=!state.cutaway;e.currentTarget.setAttribute('aria-pressed',state.cutaway);dirty=true;};
document.getElementById('reset').onclick=()=>{Object.assign(state,{yaw:.58,pitch:.30,distance:5.55});dirty=true;};
window.addEventListener('resize',()=>dirty=true);
function animate(t){if(state.orbit){state.yaw+=Math.min(t-last,40)*.00016;dirty=true;}if(dirty)redraw();last=t;requestAnimationFrame(animate);}requestAnimationFrame(animate);
}catch(e){document.getElementById('error').style.display='block';document.getElementById('error').textContent=e.message;status.textContent='Rendering error';console.error(e);}

document.getElementById('export').onclick=async()=>{
const was=state.orbit;state.orbit=false;status.textContent='Rendering high-resolution poster…';
await new Promise(r=>setTimeout(r,50));
try{const c=document.createElement('canvas');c.width=2800;c.height=1900;const ctx=c.getContext('2d');ctx.fillStyle='#080e17';ctx.fillRect(0,0,c.width,c.height);
const text=(s,x,y,size,color='#dce8e8')=>{ctx.fillStyle=color;ctx.font=`${size}px -apple-system, BlinkMacSystemFont, sans-serif`;ctx.fillText(s,x,y);};
text('VIRTUAL PLASMA BENCH  /  FIELD STUDY 01',90,85,22,'#84adae');text('Inside the field reversal.',85,188,85);text('Current and magnetic structure from an initialized FRC simulation',90,248,29,'#9aadb9');
for(let i=0;i<2;i++){let x=65+i*1370;views[i].draw(state,[1320,1180]);ctx.drawImage(views[i].canvas,x,325);text(i?'MAGNETIC FLUX DENSITY':'CURRENT DENSITY',x+35,370,31);text(i?'Poloidal field lines · color = |B|':'Toroidal current volume · color = |J|',x+35,420,22,'#9dafba');
let gradient=ctx.createLinearGradient(x+35,0,x+1280,0);for(const [f,col] of (i?[[0,'#0f40a6'],[.4,'#17e8e3'],[1,'#ffe999']]:[[0,'#470503'],[.33,'#fb3706'],[.72,'#ffbd3b'],[1,'#fff4c4']]))gradient.addColorStop(f,col);ctx.fillStyle=gradient;ctx.fillRect(x+35,1520,1245,8);
for(let k=0;k<5;k++)text((k/4*(i?DATA.bmax_mT:DATA.jmax/1000)).toFixed(2),x+35+k*296,1566,22,'#aabfc8');text(i?'mT':'kA / m²',x+1188,1610,23);}
const f=DATA.frames[state.frame];text(`SAVED TIME  ${f.time_ms} ms     INITIAL RADIUS  1 m     APPLIED BIAS  ${f.bias_mT.toFixed(3)} mT`,90,1700,26,'#b7d4d1');
text('Initialized, axisymmetric verification case — not demonstrated formation or 3D stability.',90,1756,27,'#a5b5c0');text('Fixed color scales. Illustrative opacity and field-line spacing. Heated-edge accuracy remains unresolved.',90,1802,24,'#778f9f');text('ARCHIVED CASE  source_g1_32  /  32 × 64 MATERIAL CELLS  /  EXACTLY ZERO EXTERIOR PLASMA CURRENT',90,1860,20,'#637f8d');
const blob=await new Promise(r=>c.toBlob(r,'image/png'));
if(location.protocol==='http:'||location.protocol==='https:'){let response=await fetch('/save-poster',{method:'POST',headers:{'Content-Type':'image/png'},body:blob});if(!response.ok)throw Error('Poster save failed');status.textContent='Poster saved: outputs/frc_showcase/poster.png';}
else{let a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='frc-field-study.png';a.click();status.textContent='Poster downloaded';}}
catch(e){status.textContent=e.message;}finally{state.orbit=was;dirty=true;}
};
