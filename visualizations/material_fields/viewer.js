const ids=id=>document.getElementById(id), canvases=[ids('field'),ids('current')];
let frame=DATA.initial_frame??DATA.frames.length-1,requestedFrame=frame,component=0,shared=Boolean(DATA.detailed),yaw=DATA.initial_yaw??-.65,pitch=DATA.initial_pitch??.35,zoom=DATA.initial_zoom||1,playing=false,playTimer=null,requestId=0;
let volumeRenderer=null;
if(DATA.detailed){try{volumeRenderer=createVolumeRenderer(ids('volume'));}catch(error){ids('render-warning').textContent=error.message;}}
let frameCache=new Map(),pendingFrames=new Map();
const TRAIL_US=2;let particleData=null,particleKey=null;   // optional particle trails (end of file)
function decodeVolume(v){const bytes=Uint8Array.from(atob(v.b64),c=>c.charCodeAt(0)),raw=new Float32Array(bytes.buffer),n=v.shape.reduce((a,b)=>a*b,1),k=v.channels||4;if(raw.length!==n*k)throw new Error('Incomplete current volume.');if(k===4)return raw;const out=new Float32Array(4*n),slots=v.slots||[1,3];for(let i=0;i<n;i++)for(let c=0;c<k;c++)out[4*i+slots[c]]=raw[k*i+c];return out;}
async function loadFrame(index){const cache=frameCache,pending=pendingFrames;   // a run switch replaces the maps; stale loads stay in the old ones
 const entry=DATA.frames[index];if(entry.lines){if(entry.volume&&entry.volume.b64&&!entry.volume_data)entry.volume_data=decodeVolume(entry.volume);return entry;}
 if(cache.has(index)){const value=cache.get(index);cache.delete(index);cache.set(index,value);return value;}
 if(pending.has(index))return pending.get(index);
 const task=(async()=>{const response=await fetch(entry.url,{cache:'no-store'});if(!response.ok)throw new Error('Frame could not be loaded ('+response.status+').');const value={...await response.json(),...entry};if(value.volume3d&&value.volume3d.url&&!value.volume3d_data){const raw=await fetch(value.volume3d.url,{cache:'no-store'});if(!raw.ok)throw new Error('3D volume could not be loaded ('+raw.status+').');value.volume3d_data=new Uint16Array(await raw.arrayBuffer());if(value.volume3d_data.length!==value.volume3d.shape.reduce((a,b)=>a*b,1)*4)throw new Error('Incomplete 3D volume.');}if(value.volume&&value.volume.b64&&!value.volume_data)value.volume_data=decodeVolume(value.volume);if(entry.volume){const raw=await fetch(entry.volume.url);if(!raw.ok)throw new Error('Current volume could not be loaded ('+raw.status+').');value.volume_data=new Float32Array(await raw.arrayBuffer());if(value.volume_data.length!==entry.volume.shape[0]*entry.volume.shape[1]*4)throw new Error('Incomplete current volume.');}cache.set(index,value);while(cache.size>5){const old=[...cache.keys()].find(k=>k!==frame&&k!==requestedFrame);cache.delete(old);}return value;})();
 pending.set(index,task);try{return await task;}finally{pending.delete(index);}
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
function color(t,signed=false){t=Math.max(0,Math.min(1,t));if(DATA.detailed&&!signed){const stops=[[15,31,140],[0,179,217],[64,217,97],[255,224,38],[255,77,15],[199,6,4]],v=t*5,i=Math.min(4,Math.floor(v));return rgb(mix(stops[i],stops[i+1],v-i));}return rgb(t<.5?mix(signed?[49,90,158]:[52,120,186],signed?[241,245,245]:[103,221,206],2*t):mix(signed?[241,245,245]:[103,221,206],signed?[211,90,71]:[255,241,173],2*t-1));}
function project(p,c){const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);const x=cy*p[0]+sy*p[1],y=-sy*p[0]+cy*p[1],z=cp*p[2]-sp*y,depth=sp*p[2]+cp*y;const f=zoom*Math.min(c.width,c.height)*.37/(DATA.detailed?1:1+depth*.18);return [c.width/2+x*f,c.height/2-z*f,depth];}
function line(ctx,points,c,stroke,width=1){if(points.length<2)return;ctx.strokeStyle=stroke;ctx.lineWidth=width;ctx.beginPath();points.forEach((p,i)=>{const s=project(p,c);if(i)ctx.lineTo(s[0],s[1]);else ctx.moveTo(s[0],s[1]);});ctx.stroke();}
function draw(){const f=frameCache.get(frame)||DATA.frames[frame];if(!f.lines)return;const bmax=shared?DATA.bmax:f.bmax,jmax=shared?DATA.jmax[component]:f.jmax[component];
 canvases.forEach((c,index)=>{const ratio=Math.min(window.devicePixelRatio||1,DATA.detailed?1.5:Infinity),w=Math.round(c.clientWidth*ratio),h=Math.round(c.clientHeight*ratio);if(c.width!==w||c.height!==h){c.width=w;c.height=h;}const ctx=c.getContext('2d');ctx.clearRect(0,0,w,h);if(!(DATA.detailed&&index===1&&volumeRenderer)){ctx.fillStyle=DATA.detailed?'#000':'#0e202d';ctx.fillRect(0,0,w,h);}
  if(!DATA.detailed||ids('geometry').checked){f.outline.forEach(p=>line(ctx,p,c,DATA.detailed?'#253138':'#527080',.7*ratio));
  (DATA.chamber||[]).forEach(p=>line(ctx,p,c,DATA.detailed?'#36414a':'#7c8a94',.6*ratio));
  if(DATA.coil_groups)DATA.coil_groups.forEach(g=>g.paths.forEach(p=>line(ctx,p,c,g.color,1.8*ratio)));
  else (DATA.coils||[]).forEach(p=>line(ctx,p,c,'#db9160',1.8*ratio));}
  if(index===0){const strokes=new Map(),copies=DATA.detailed&&!DATA.no_copies?Number(ids('density').value):1;for(let copy=0;copy<copies;copy++){const angle=copy*Math.PI/(2*copies),ca=Math.cos(angle),sa=Math.sin(angle);f.lines.forEach(path=>{const projected=path.points.map(p=>project(copy?[p[0],p[1]*ca-p[2]*sa,p[1]*sa+p[2]*ca]:p,c));for(let i=1;i<projected.length;i++){const level=Math.round(255*Math.max(0,Math.min(1,bmax>0?path.field[i]/bmax:0)));if(!strokes.has(level))strokes.set(level,new Path2D());const stroke=strokes.get(level);stroke.moveTo(projected[i-1][0],projected[i-1][1]);stroke.lineTo(projected[i][0],projected[i][1]);}});}ctx.lineWidth=(DATA.detailed?.65:1.2)*ratio;ctx.globalAlpha=DATA.detailed?.65:1;for(const [level,stroke] of strokes){ctx.strokeStyle=color(level/255);ctx.stroke(stroke);}ctx.globalAlpha=1;if(particleData&&ids('particles')&&ids('particles').checked)drawParticles(ctx,c,ratio,f.time_us);}
  else if(volumeRenderer&&(f.volume_data||f.volume3d_data)){const vc=ids('volume');if(vc.width!==w||vc.height!==h){vc.width=w;vc.height=h;}volumeRenderer.draw(f,{yaw,pitch,zoom,maximum:jmax,component,opacity:Number(ids('opacity').value),threshold:Number(ids('threshold').value)});}
  else{const points=f.current.map(p=>({p:project(p.slice(0,3),c),j:p[3+component]})).sort((a,b)=>b.p[2]-a.p[2]);ctx.globalAlpha=.52;for(const p of points){const scaled=jmax>0?p.j/jmax:0,t=component?.5+.5*scaled:scaled;ctx.fillStyle=color(t,component!==0);ctx.beginPath();ctx.arc(p.p[0],p.p[1],1.9*ratio,0,2*Math.PI);ctx.fill();}ctx.globalAlpha=1;}
  if(!ids('separatrix')||ids('separatrix').checked)(f.separatrix||[]).forEach(p=>line(ctx,p,c,DATA.separatrix_color||'#f2c65c',(index?1.1:1.6)*ratio));
  if(f.ions&&ids('ions')&&ids('ions').checked){const sp=DATA.ion_speed_max||1;ctx.globalAlpha=.75;for(const q of f.ions){const s=project(q,c);ctx.fillStyle=color(Math.min(1,q[3]/sp));ctx.fillRect(s[0]-ratio,s[1]-ratio,2*ratio,2*ratio);}ctx.globalAlpha=1;}
  if(DATA.detailed){const length=.05,bar=length/DATA.physical_scale_m*zoom*Math.min(w,h)*.37;ctx.strokeStyle='#788b97';ctx.lineWidth=ratio;ctx.beginPath();ctx.moveTo(23*ratio,h-23*ratio);ctx.lineTo(23*ratio+bar,h-23*ratio);ctx.stroke();ctx.fillStyle='#8c9da9';ctx.font=10*ratio+'px system-ui';ctx.fillText('50 mm',23*ratio,h-31*ratio);}
 });
 const timeLabel=DATA.display_interval_ns?(f.time_us*1000).toFixed(2)+' ns':f.time_us.toFixed(3)+' µs';ids('time').value=frame;ids('time-label').textContent=timeLabel;ids('metric-time').textContent=timeLabel;ids('field-max').textContent=(bmax*(DATA.detailed?1e3:1e6)).toPrecision(3)+(DATA.detailed?' mT':' µT');
 const units=DATA.detailed&&jmax>=1e6?'MA/m²':jmax>=1000?'kA/m²':'A/m²',scale=DATA.detailed&&jmax>=1e6?1e6:jmax>=1000?1000:1;ids('current-min').textContent=(component?-jmax/scale:0).toPrecision(3)+' '+units;ids('current-max').textContent=(jmax/scale).toPrecision(3)+' '+units;
 ids('current-title').textContent=['Current density magnitude','Radial current density','Azimuthal current density','Axial current density'][component];ids('current-gradient').style.background=component?(DATA.detailed?'linear-gradient(90deg,#265ff2,#f0f0f0,#fa3017)':'linear-gradient(90deg,#315a9e,#f1f5f5,#d35a47)'):DATA.detailed?'linear-gradient(90deg,#0f1f8c,#00b3d9,#40d961,#ffe026,#ff4d0f,#c70604)':'linear-gradient(90deg,#3478ba,#67ddce,#fff1ad)';
 ids('metric-field').textContent=(f.bmax*(DATA.detailed?1e3:1e6)).toPrecision(4)+(DATA.detailed?' mT':' µT');ids('metric-current').textContent=DATA.detailed?(f.jmax[0]/(f.jmax[0]>=1e6?1e6:f.jmax[0]>=1000?1000:1)).toPrecision(4)+(f.jmax[0]>=1e6?' MA/m²':f.jmax[0]>=1000?' kA/m²':' A/m²'):f.jmax[0].toPrecision(4)+' A/m²';ids('metric-temp').textContent=(f.temperature||[]).map(v=>v==null?'—':DATA.temperature_digits?v.toPrecision(DATA.temperature_digits):v.toFixed(3)).join(' / ')+' eV';ids('scale').textContent=shared?'Use per-state color scale':'Use shared color scale';ids('status').textContent='Saved state '+(frame+1)+' of '+DATA.frames.length+' · '+(shared?'shared':'per-state')+' color scale · '+(DATA.status_note||'initialized flux, no formation claim');
 if(ids('timing')){const dt=DATA.solver_dt_s,gap=frame?(f.time_us-DATA.frames[frame-1].time_us)*1000:0,completed=DATA.solver_completed_steps??DATA.solver_steps;const steps=completed===DATA.solver_steps?DATA.solver_steps?.toLocaleString()+' steps':completed.toLocaleString()+' / '+DATA.solver_steps.toLocaleString()+' completed steps';const sampling=DATA.display_interval_ns?' · '+DATA.display_interval_ns+' ns target, nearest saved states':'';ids('timing').textContent=dt===undefined?'Display: '+DATA.frames.length+' selected saved states.':'Simulation step: '+(dt*1e9).toFixed(3)+' ns · '+steps+'. Display: '+DATA.frames.length+' states'+sampling+(frame?' · gap from preceding frame '+gap.toFixed(2)+' ns.':' · initial state.');}
 if(ids('drive'))ids('drive').textContent=f.drive_center_t===undefined?'':(DATA.drive_label||'Prescribed vacuum-centre field')+': '+(f.drive_center_t*1e3).toFixed(3)+' mT';
 if(ids('drive-groups')){const el=ids('drive-groups');el.replaceChildren();(f.drive_groups||[]).forEach((g,i)=>{const span=document.createElement('span');span.style.color=DATA.coil_groups[i].color;span.textContent=(i?' · ':'')+g.name+': '+(g.center_t*1e3).toFixed(3)+' mT';el.append(span);});}
 if(ids('frame-summary'))ids('frame-summary').textContent=f.summary||'';
 if(ids('ion-legend'))ids('ion-legend').textContent=f.ions?.length&&ids('ions')?.checked?'Ion dot colors: speed, 0–'+((DATA.ion_speed_max||1)/1000).toPrecision(3)+' km/s; faster ions saturate. Each frame is a new sample.':'';
 if(ids('motion'))ids('motion').textContent=f.midplane_radius_m===undefined?'':(DATA.motion_label||'Plasma radius — middle')+': '+(f.midplane_radius_m*1000).toFixed(2)+' mm · maximum: '+(f.maximum_radius_m*1000).toFixed(2)+' mm';}
let scrubTimer=null;ids('time').max=DATA.frames.length-1;ids('time').value=frame;ids('time').addEventListener('input',e=>{stopPlayback();requestId++;requestedFrame=+e.target.value;clearTimeout(scrubTimer);scrubTimer=setTimeout(()=>showFrame(requestedFrame),65);});ids('component').addEventListener('change',e=>{component=+e.target.value;draw();});ids('scale').onclick=()=>{shared=!shared;draw();};ids('reset').onclick=()=>{yaw=DATA.initial_yaw??-.65;pitch=DATA.initial_pitch??.35;zoom=DATA.initial_zoom||1;if(ids('geometry'))ids('geometry').checked=false;draw();};
for(const id of ['opacity','threshold','density','geometry','separatrix','ions'])if(ids(id))ids(id).addEventListener('input',draw);
if(ids('apparatus'))ids('apparatus').onclick=()=>{yaw=-.65;pitch=.35;zoom=1.65;ids('geometry').checked=true;draw();};
if(ids('previous'))ids('previous').onclick=()=>stepFrame(-1);if(ids('next'))ids('next').onclick=()=>stepFrame(1);if(ids('speed'))ids('speed').value=DATA.playback_fps||1;
ids('play').onclick=()=>{if(playing){stopPlayback();requestId++;requestedFrame=frame;draw();}else{clearTimeout(scrubTimer);playing=true;ids('play').textContent='Pause';advance();}};
canvases.forEach(c=>{let drag=null;c.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY];c.setPointerCapture(e.pointerId);});c.addEventListener('pointermove',e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.008;pitch=Math.max(-1.45,Math.min(1.45,pitch+(e.clientY-drag[1])*.008));drag=[e.clientX,e.clientY];draw();});c.addEventListener('pointerup',()=>drag=null);c.addEventListener('pointercancel',()=>drag=null);c.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.4,Math.min(DATA.detailed?10:3,zoom*Math.exp(-e.deltaY*.001)));draw();},{passive:false});});
ids('scope').textContent=DATA.scope;window.addEventListener('resize',draw);showFrame(frame);
function showRunText(){if(ids('run-title'))ids('run-title').textContent=DATA.title||'';if(ids('run-description'))ids('run-description').textContent=DATA.description||'';if(DATA.title)document.title=DATA.title+' · FRC formation';ids('scope').textContent=DATA.scope||'';
 if(ids('model-badge'))ids('model-badge').textContent=DATA.model_label||'Axisymmetric · 3D rendering';
 if(ids('temperature-label'))ids('temperature-label').textContent=DATA.temperature_label||'Te / Ti inside separatrix (Te = Ti)';
 if(ids('model-notice'))ids('model-notice').textContent=DATA.model_notice||'Saved axisymmetric formation states. Three-dimensional stability and confinement are not established.';
 if(ids('render-method'))ids('render-method').textContent=DATA.rendering_note||'The volume is revolved from the saved axisymmetric current. Colors show A/m²; opacity is illustrative, not emitted light. Field lines are interpolated paths, not particle trajectories. All directions share one length scale.';
 if(ids('density')){ids('density').hidden=!!DATA.no_copies;const label=document.querySelector('label[for="density"]');if(label)label.hidden=!!DATA.no_copies;}
 if(ids('ions'))ids('ions').disabled=!DATA.no_copies;
 for(const [i,label] of [[1,'Radial Jr'],[3,'Axial Jz']]){const option=ids('component').querySelector('option[value="'+i+'"]');if(option){option.disabled=!(DATA.jmax[i]>0);option.textContent=label+(option.disabled?' (zero in this run)':'');}}
 if(ids('component').selectedOptions[0]?.disabled){component=0;ids('component').value=0;}}
async function switchRun(name){
 const entry=(DATA.runs||[]).find(r=>r.name===name);if(!entry||name===DATA.run)return;
 const wasPlaying=playing;stopPlayback();clearTimeout(scrubTimer);const request=++requestId,t=DATA.frames[frame]?DATA.frames[frame].time_us:0;
 ids('status').textContent='Loading '+entry.label+'…';let next;
 try{const response=await fetch(entry.manifest,{cache:'no-store'});if(!response.ok)throw new Error('Run could not be loaded ('+response.status+').');next=await response.json();}
 catch(error){if(request===requestId){ids('status').textContent='Unable to load this run. '+error.message;if(ids('run'))ids('run').value=DATA.run;}return;}
 if(request!==requestId)return;
 const runs=DATA.runs;for(const k of Object.keys(DATA))delete DATA[k];Object.assign(DATA,next,{runs,run:name});frameCache=new Map();pendingFrames=new Map();
 let best=0;DATA.frames.forEach((f,i)=>{if(Math.abs(f.time_us-t)<Math.abs(DATA.frames[best].time_us-t))best=i;});   // keep the simulated time; camera untouched
 ids('time').max=DATA.frames.length-1;frame=requestedFrame=best;showRunText();particleData=null;particleKey=null;setupParticles();if(ids('run'))ids('run').value=name;
 try{history.replaceState(null,'','?run='+encodeURIComponent(name));}catch(e){}
 await showFrame(best);if(wasPlaying&&request+1===requestId){playing=true;ids('play').textContent='Pause';advance();}
}
if(DATA.runs&&ids('run')){const sel=ids('run');DATA.runs.forEach(r=>{const o=document.createElement('option');o.value=r.name;o.textContent=r.label;sel.append(o);});
 sel.value=DATA.run;sel.addEventListener('change',e=>switchRun(e.target.value));showRunText();
 const wanted=new URLSearchParams(location.search).get('run');if(wanted&&wanted!==DATA.run)switchRun(wanted);}
// Optional particle trails (examples/particle_tracking.py via formation_display): DATA.particles is the url of
// actual integrated orbits, sampled every sample_us from t0_us; each frame shows the last TRAIL_US of each orbit.
function particleSet(){const d=particleData;return d?d.sets[Math.min(+(ids('particle-set')?.value||0),d.sets.length-1)]:null;}
function particleLegend(text){const el=ids('particle-legend');if(!el)return;el.replaceChildren();if(text!==undefined){el.textContent=text;return;}const d=particleData;if(!d||!ids('particles').checked)return;
 const head=document.createElement('span');head.textContent=d.species+' orbits ('+particleSet().label+'), last '+TRAIL_US+' µs of each, from '+d.t0_us+' µs: ';el.append(head);
 [...new Set(particleSet().trails.map(t=>t.class))].forEach(k=>{const s=document.createElement('span');s.style.color=d.colors[k]||'#fff';s.textContent='● '+k+'  ';el.append(s);});
 const x=document.createElement('span');x.style.color='#ff5a4f';x.textContent='× lost at wall or end';el.append(x);}
async function loadParticles(){const url=DATA.particles;if(!url||particleKey===url)return;particleKey=url;particleLegend('Loading particle orbits…');
 try{const r=await fetch(url,{cache:'no-store'});if(!r.ok)throw new Error('('+r.status+')');const d=await r.json();if(particleKey!==url)return;particleData=d;
  const sel=ids('particle-set');if(sel){sel.replaceChildren();d.sets.forEach((s,i)=>{const o=document.createElement('option');o.value=i;o.textContent=s.label;sel.append(o);});sel.hidden=d.sets.length<2;}
  particleLegend();draw();}
 catch(error){if(particleKey===url){particleKey=null;particleLegend('Particle orbits could not be loaded '+error.message);}}}
function drawParticles(ctx,c,ratio,t){const d=particleData,set=particleSet(),span=Math.round(TRAIL_US/d.sample_us);
 for(const tr of set.trails){const n=tr.points.length,lost=tr.lost_us!=null&&t>=tr.lost_us;let end=lost?n-1:Math.min(tr.lost_us!=null?n-2:n-1,Math.floor((t-d.t0_us)/d.sample_us+1e-6));
  if(end<1)continue;const start=Math.max(0,end-span),pts=tr.points.slice(start,end+1).map(p=>project(p,c));
  ctx.strokeStyle=d.colors[tr.class]||'#fff';ctx.globalAlpha=lost?.45:.95;ctx.lineWidth=1.3*ratio;ctx.beginPath();pts.forEach((s,i)=>i?ctx.lineTo(s[0],s[1]):ctx.moveTo(s[0],s[1]));ctx.stroke();
  const h=pts[pts.length-1];ctx.globalAlpha=1;if(lost){ctx.strokeStyle='#ff5a4f';ctx.lineWidth=1.6*ratio;const a=4*ratio;ctx.beginPath();ctx.moveTo(h[0]-a,h[1]-a);ctx.lineTo(h[0]+a,h[1]+a);ctx.moveTo(h[0]+a,h[1]-a);ctx.lineTo(h[0]-a,h[1]+a);ctx.stroke();}
  else{ctx.fillStyle=d.colors[tr.class]||'#fff';ctx.beginPath();ctx.arc(h[0],h[1],2.2*ratio,0,2*Math.PI);ctx.fill();}}
 ctx.globalAlpha=1;}
function setupParticles(){const box=ids('particles');if(!box)return;box.disabled=!DATA.particles;if(ids('particle-set'))ids('particle-set').hidden=true;
 if(!DATA.particles){particleLegend(box.checked?'No tracked orbits for this run.':'');return;}if(box.checked)loadParticles();else particleLegend('');}
if(ids('particles')){ids('particles').addEventListener('input',()=>{if(ids('particles').checked&&!particleData)setupParticles();particleLegend();draw();});
 if(ids('particle-set'))ids('particle-set').addEventListener('input',()=>{particleLegend();draw();});setupParticles();}
