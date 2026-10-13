const state={yaw:.50,pitch:.30,distance:3.70,density:1.05,cutaway:true,frame:DATA.frames.length-1,orbit:false,localScale:false};
let views=[],dirty=true,last=0;
const status=document.getElementById('status');
const number=(x,d=2)=>x.toLocaleString('en-US',{minimumFractionDigits:d,maximumFractionDigits:d});
function currentMax(){return state.localScale?DATA.frames[state.frame].sampled_Jmax_A_m2:DATA.jmax;}
function metrics(){const f=DATA.frames[state.frame];
document.getElementById('time-readout').textContent=`${number(f.time_us)} µs / ${f.step}`;
document.getElementById('current-readout').textContent=`${number(f.sampled_Jmax_A_m2/1000,3)} kA/m²`;
document.getElementById('magnetic-readout').textContent=`${number(f.sampled_Bmin_mT,4)}–${number(f.sampled_Bmax_mT,4)} mT`;
document.getElementById('axial-readout').textContent=`+${number(f.sampled_Bzmin_mT,4)} mT`;
for(const [id,min,max,digits] of [['j-ticks',0,currentMax()/1000,state.localScale&&currentMax()<1000?4:1],['b-ticks',DATA.bmin_mT,DATA.bmax_mT,3]])document.getElementById(id).innerHTML=[0,.25,.5,.75,1].map(t=>`<span>${number(min+t*(max-min),digits)}</span>`).join('');
document.getElementById('current-scale-label').textContent=state.localScale?'Per-frame current scale · changes with time':'Shared current scale across saved times';
document.getElementById('scale-note').textContent=state.localScale?'Per-frame scale reveals weak structure; compare values using the legends.':'Shared scale preserves the change in current magnitude.';}
function selectFrame(i){state.frame=i;for(const v of views)v.setFrame(i);document.querySelectorAll('[data-frame]').forEach(x=>{const selected=+x.dataset.frame===i;x.classList.toggle('selected',selected);x.setAttribute('aria-pressed',selected);});metrics();dirty=true;}
try{
document.getElementById('times').innerHTML=DATA.frames.map((f,i)=>`<button data-frame="${i}" aria-pressed="${i===state.frame}" class="${i===state.frame?'selected':''}">${number(f.time_us)} µs${i===DATA.frames.length-1?' · latest':''}</button>`).join('');
views=[new FieldView(document.getElementById('current'),DATA,0),new FieldView(document.getElementById('magnetic'),DATA,1)];metrics();
for(const canvas of document.querySelectorAll('canvas')){let drag=null;
canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);canvas.style.cursor='grabbing';};
canvas.onpointermove=e=>{if(!drag)return;state.yaw-=(e.clientX-drag[0])*.006;state.pitch=Math.max(-1.25,Math.min(1.25,state.pitch+(e.clientY-drag[1])*.006));drag=[e.clientX,e.clientY];dirty=true;};
const stop=()=>{drag=null;canvas.style.cursor='grab';};canvas.onpointerup=stop;canvas.onpointercancel=stop;
canvas.addEventListener('wheel',e=>{e.preventDefault();state.distance=Math.max(2.8,Math.min(7,state.distance*Math.exp(e.deltaY*.001)));dirty=true;},{passive:false});}
document.querySelectorAll('[data-frame]').forEach(b=>b.onclick=()=>selectFrame(+b.dataset.frame));
document.getElementById('orbit').onclick=e=>{state.orbit=!state.orbit;e.currentTarget.setAttribute('aria-pressed',state.orbit);};
document.getElementById('cutaway').onclick=e=>{state.cutaway=!state.cutaway;e.currentTarget.setAttribute('aria-pressed',state.cutaway);dirty=true;};
document.getElementById('scale').onclick=e=>{state.localScale=!state.localScale;e.currentTarget.setAttribute('aria-pressed',state.localScale);metrics();dirty=true;};
document.getElementById('opacity').oninput=e=>{state.density=+e.target.value;dirty=true;};
document.getElementById('reset').onclick=()=>{Object.assign(state,{yaw:.50,pitch:.30,distance:3.70,orbit:false});document.getElementById('orbit').setAttribute('aria-pressed','false');dirty=true;};
window.addEventListener('resize',()=>dirty=true);
function animate(t){if(state.orbit){state.yaw+=Math.min(t-last,40)*.00016;dirty=true;}if(dirty){for(const v of views)v.draw(state);dirty=false;}last=t;requestAnimationFrame(animate);}requestAnimationFrame(animate);
}catch(e){document.getElementById('error').style.display='block';document.getElementById('error').textContent=e.message;status.textContent='Rendering error';console.error(e);}

async function poster(){const orbit=state.orbit;state.orbit=false;status.textContent='Rendering PNG…';
try{const c=document.createElement('canvas');c.width=2800;c.height=1900;const ctx=c.getContext('2d');ctx.fillStyle='#08101a';ctx.fillRect(0,0,c.width,c.height);
const text=(s,x,y,size,color='#dce8e8')=>{ctx.fillStyle=color;ctx.font=`${size}px -apple-system,BlinkMacSystemFont,sans-serif`;ctx.fillText(s,x,y);};
const f=DATA.frames[state.frame];
text('VIRTUAL PLASMA BENCH / CURRENT SIMULATION',90,90,24,'#82adae');text('The saved magnetic and current fields.',85,195,76);text(`Saved time ${number(f.time_us)} µs · step ${f.step} · 70 mm diameter × 240 mm · 5 mT axial bias`,90,263,30,'#9fbac6');
for(let i=0;i<2;i++){const x=65+i*1370;views[i].draw(state,[1320,1100]);ctx.drawImage(views[i].canvas,x,400);text(i?'MAGNETIC FIELD LINES':'CURRENT DENSITY',x+30,350,31);text(i?'Full vector B · actual deflection · color = |B|':'All components · volume color = |J| · illustrative opacity',x+30,393,23,'#9db4c1');
let gradient=ctx.createLinearGradient(x+35,0,x+1280,0);for(const [t,col] of (i?[[0,'#0f40a6'],[.4,'#17e8e3'],[1,'#ffe999']]:[[0,'#470503'],[.4,'#ff2e06'],[.78,'#ffb83d'],[1,'#fff4c4']]))gradient.addColorStop(t,col);ctx.fillStyle=gradient;ctx.fillRect(x+35,1510,1245,8);
const lo=i?DATA.bmin_mT:0,hi=i?DATA.bmax_mT:currentMax()/1000;
for(let k=0;k<5;k++)text(number(lo+(hi-lo)*k/4,i?3:(hi<1?4:1)),x+35+k*287,1560,24,'#acbfca');text(i?'mT':`kA/m² · ${state.localScale?'per-frame':'shared'} scale`,x+35,1612,23,'#8daab6');}
text('Run stopped at 28.33 µs: step 18 failed to converge. Requested duration: 100 µs.',90,1694,32,'#e4c28f');
text('Contact verification fixture · positive axial field · no FRC formation demonstrated.',90,1754,30,'#bdced8');
text('Axisymmetric data revolved in 3D · periodic ends · line spacing and opacity are presentation choices.',90,1812,25,'#8ba6b4');
text('SOURCE 6a96344a6400 / 17 ACCEPTED STEPS AND REJECTION REPRODUCED IN COMPLETE REPLAY',90,1863,21,'#6f8d9e');
const blob=await new Promise(resolve=>c.toBlob(resolve,'image/png'));if(!blob)throw Error('PNG encoding failed');
if(location.protocol==='http:'){const response=await fetch('/save-poster',{method:'POST',headers:{'Content-Type':'image/png'},body:blob});if(!response.ok)throw Error('PNG save failed');status.textContent='Saved poster.png in outputs/trajectory_display';}
else{const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='frc-current-simulation.png';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);status.textContent='PNG downloaded';}
}catch(e){status.textContent=e.message;console.error(e);}finally{state.orbit=orbit;dirty=true;}}
document.getElementById('export').onclick=poster;
