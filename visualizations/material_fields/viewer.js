const ids=id=>document.getElementById(id), canvases=[ids('field'),ids('current')];
let frame=DATA.initial_frame??DATA.frames.length-1,requestedFrame=frame,component=0,shared=false,yaw=-.65,pitch=.35,zoom=DATA.initial_zoom||1,playing=false,playTimer=null,requestId=0;
const frameCache=new Map(),pendingFrames=new Map();
async function loadFrame(index){
 const entry=DATA.frames[index];if(entry.lines)return entry;
 if(frameCache.has(index)){const value=frameCache.get(index);frameCache.delete(index);frameCache.set(index,value);return value;}
 if(pendingFrames.has(index))return pendingFrames.get(index);
 const task=(async()=>{const response=await fetch(entry.url,{cache:'no-store'});if(!response.ok)throw new Error('Frame could not be loaded ('+response.status+').');const value=await response.json();frameCache.set(index,value);while(frameCache.size>5){const old=[...frameCache.keys()].find(k=>k!==frame&&k!==requestedFrame);frameCache.delete(old);}return value;})();
 pendingFrames.set(index,task);try{return await task;}finally{pendingFrames.delete(index);}
}
async function showFrame(index){
 requestedFrame=Math.max(0,Math.min(DATA.frames.length-1,index));const target=requestedFrame,request=++requestId;
 ids('time').value=target;ids('status').textContent='Loading saved state '+(target+1)+' of '+DATA.frames.length+'…';
 try{await loadFrame(target);if(request!==requestId)return;frame=target;draw();if(target+1<DATA.frames.length)loadFrame(target+1).catch(()=>{});}
 catch(error){if(request!==requestId)return;stopPlayback();ids('status').textContent='Unable to load this saved state. '+error.message;ids('time').value=frame;}
}
function stopPlayback(){playing=false;clearTimeout(playTimer);ids('play').textContent='Play saved states';}
async function advance(){
 if(!playing)return;const start=performance.now();await showFrame((frame+1)%DATA.frames.length);
 if(playing)playTimer=setTimeout(advance,Math.max(0,1000/Number(ids('speed')?.value||DATA.playback_fps||1)-(performance.now()-start)));
}
function stepFrame(delta){clearTimeout(scrubTimer);stopPlayback();showFrame(requestedFrame+delta);}
const mix=(a,b,t)=>a.map((v,i)=>Math.round(v+(b[i]-v)*t));
const rgb=c=>'rgb('+c.join(',')+')';
function color(t,signed=false){t=Math.max(0,Math.min(1,t));return rgb(t<.5?mix(signed?[49,90,158]:[52,120,186],signed?[241,245,245]:[103,221,206],2*t):mix(signed?[241,245,245]:[103,221,206],signed?[211,90,71]:[255,241,173],2*t-1));}
function project(p,c){const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);const x=cy*p[0]+sy*p[1],y=-sy*p[0]+cy*p[1],z=cp*p[2]-sp*y,depth=sp*p[2]+cp*y;const f=zoom*Math.min(c.width,c.height)*.37/(1+depth*.18);return [c.width/2+x*f,c.height/2-z*f,depth];}
function line(ctx,points,c,stroke,width=1){if(points.length<2)return;ctx.strokeStyle=stroke;ctx.lineWidth=width;ctx.beginPath();points.forEach((p,i)=>{const s=project(p,c);if(i)ctx.lineTo(s[0],s[1]);else ctx.moveTo(s[0],s[1]);});ctx.stroke();}
function draw(){const f=frameCache.get(frame)||DATA.frames[frame];if(!f.lines)return;const bmax=shared?DATA.bmax:f.bmax,jmax=shared?DATA.jmax[component]:f.jmax[component];
 canvases.forEach((c,index)=>{const ratio=window.devicePixelRatio||1,w=Math.round(c.clientWidth*ratio),h=Math.round(c.clientHeight*ratio);if(c.width!==w||c.height!==h){c.width=w;c.height=h;}const ctx=c.getContext('2d');ctx.clearRect(0,0,w,h);ctx.fillStyle='#0e202d';ctx.fillRect(0,0,w,h);
  f.outline.forEach(p=>line(ctx,p,c,'#527080',.7*ratio));
  (DATA.chamber||[]).forEach(p=>line(ctx,p,c,'#7c8a94',.6*ratio));
  if(DATA.coil_groups)DATA.coil_groups.forEach(g=>g.paths.forEach(p=>line(ctx,p,c,g.color,1.8*ratio)));
  else (DATA.coils||[]).forEach(p=>line(ctx,p,c,'#db9160',1.8*ratio));
  if(index===0){const strokes=new Map();f.lines.forEach(path=>{const projected=path.points.map(p=>project(p,c));for(let i=1;i<projected.length;i++){const level=Math.round(255*Math.max(0,Math.min(1,bmax>0?path.field[i]/bmax:0)));if(!strokes.has(level))strokes.set(level,new Path2D());const stroke=strokes.get(level);stroke.moveTo(projected[i-1][0],projected[i-1][1]);stroke.lineTo(projected[i][0],projected[i][1]);}});ctx.lineWidth=1.2*ratio;for(const [level,stroke] of strokes){ctx.strokeStyle=color(level/255);ctx.stroke(stroke);}}
  else{const points=f.current.map(p=>({p:project(p.slice(0,3),c),j:p[3+component]})).sort((a,b)=>b.p[2]-a.p[2]);ctx.globalAlpha=.52;for(const p of points){const scaled=jmax>0?p.j/jmax:0,t=component?.5+.5*scaled:scaled;ctx.fillStyle=color(t,component!==0);ctx.beginPath();ctx.arc(p.p[0],p.p[1],1.9*ratio,0,2*Math.PI);ctx.fill();}ctx.globalAlpha=1;}
 });
 const timeLabel=DATA.display_interval_ns?(f.time_us*1000).toFixed(2)+' ns':f.time_us.toFixed(3)+' µs';ids('time').value=frame;ids('time-label').textContent=timeLabel;ids('metric-time').textContent=timeLabel;ids('field-max').textContent=(bmax*1e6).toPrecision(3)+' µT';
 const units=jmax>=1000?'kA/m²':'A/m²',scale=jmax>=1000?1000:1;ids('current-min').textContent=(component?-jmax/scale:0).toPrecision(3)+' '+units;ids('current-max').textContent=(jmax/scale).toPrecision(3)+' '+units;
 ids('current-title').textContent=['Current density magnitude','Radial current density','Azimuthal current density','Axial current density'][component];ids('current-gradient').style.background=component?'linear-gradient(90deg,#315a9e,#f1f5f5,#d35a47)':'linear-gradient(90deg,#3478ba,#67ddce,#fff1ad)';
 ids('metric-field').textContent=(f.bmax*1e6).toPrecision(4)+' µT';ids('metric-current').textContent=f.jmax[0].toPrecision(4)+' A/m²';ids('metric-temp').textContent=f.temperature.map(v=>v.toFixed(3)).join(' / ')+' eV';ids('scale').textContent=shared?'Use per-state color scale':'Use shared color scale';ids('status').textContent='Saved state '+(frame+1)+' of '+DATA.frames.length+' · '+(shared?'shared':'per-state')+' color scale · '+(DATA.status_note||'initialized flux, no formation claim');
 if(ids('timing')){const dt=DATA.solver_dt_s,gap=frame?(f.time_us-DATA.frames[frame-1].time_us)*1000:0,completed=DATA.solver_completed_steps??DATA.solver_steps;const steps=completed===DATA.solver_steps?DATA.solver_steps?.toLocaleString()+' steps':completed.toLocaleString()+' / '+DATA.solver_steps.toLocaleString()+' completed steps';const sampling=DATA.display_interval_ns?' · '+DATA.display_interval_ns+' ns target, nearest saved states':'';ids('timing').textContent=dt===undefined?'Display: '+DATA.frames.length+' selected saved states.':'Simulation step: '+(dt*1e9).toFixed(3)+' ns · '+steps+'. Display: '+DATA.frames.length+' states'+sampling+(frame?' · gap from preceding frame '+gap.toFixed(2)+' ns.':' · initial state.');}
 if(ids('drive'))ids('drive').textContent=f.drive_center_t===undefined?'':'Prescribed vacuum-centre field: '+(f.drive_center_t*1e3).toFixed(3)+' mT';
 if(ids('drive-groups')){const el=ids('drive-groups');el.replaceChildren();(f.drive_groups||[]).forEach((g,i)=>{const span=document.createElement('span');span.style.color=DATA.coil_groups[i].color;span.textContent=(i?' · ':'')+g.name+': '+(g.center_t*1e3).toFixed(3)+' mT';el.append(span);});}
 if(ids('motion'))ids('motion').textContent=f.midplane_radius_m===undefined?'':'Plasma radius — middle: '+(f.midplane_radius_m*1000).toFixed(2)+' mm · maximum: '+(f.maximum_radius_m*1000).toFixed(2)+' mm';}
let scrubTimer=null;ids('time').max=DATA.frames.length-1;ids('time').value=frame;ids('time').addEventListener('input',e=>{stopPlayback();requestId++;requestedFrame=+e.target.value;clearTimeout(scrubTimer);scrubTimer=setTimeout(()=>showFrame(requestedFrame),65);});ids('component').addEventListener('change',e=>{component=+e.target.value;draw();});ids('scale').onclick=()=>{shared=!shared;draw();};ids('reset').onclick=()=>{yaw=-.65;pitch=.35;zoom=DATA.initial_zoom||1;draw();};
if(ids('previous'))ids('previous').onclick=()=>stepFrame(-1);if(ids('next'))ids('next').onclick=()=>stepFrame(1);if(ids('speed'))ids('speed').value=DATA.playback_fps||1;
ids('play').onclick=()=>{if(playing){stopPlayback();requestId++;requestedFrame=frame;draw();}else{clearTimeout(scrubTimer);playing=true;ids('play').textContent='Pause';advance();}};
canvases.forEach(c=>{let drag=null;c.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY];c.setPointerCapture(e.pointerId);});c.addEventListener('pointermove',e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.008;pitch=Math.max(-1.45,Math.min(1.45,pitch+(e.clientY-drag[1])*.008));drag=[e.clientX,e.clientY];draw();});c.addEventListener('pointerup',()=>drag=null);c.addEventListener('pointercancel',()=>drag=null);c.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.4,Math.min(3,zoom*Math.exp(-e.deltaY*.001)));draw();},{passive:false});});
ids('scope').textContent=DATA.scope;window.addEventListener('resize',draw);showFrame(frame);
