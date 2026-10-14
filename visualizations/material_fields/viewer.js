const ids=id=>document.getElementById(id), canvases=[ids('field'),ids('current')];
let frame=DATA.frames.length-1,component=0,shared=false,yaw=-.65,pitch=.35,zoom=DATA.initial_zoom||1,playing=null;
const mix=(a,b,t)=>a.map((v,i)=>Math.round(v+(b[i]-v)*t));
const rgb=c=>'rgb('+c.join(',')+')';
function color(t,signed=false){t=Math.max(0,Math.min(1,t));return rgb(t<.5?mix(signed?[49,90,158]:[52,120,186],signed?[241,245,245]:[103,221,206],2*t):mix(signed?[241,245,245]:[103,221,206],signed?[211,90,71]:[255,241,173],2*t-1));}
function project(p,c){const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);const x=cy*p[0]+sy*p[1],y=-sy*p[0]+cy*p[1],z=cp*p[2]-sp*y,depth=sp*p[2]+cp*y;const f=zoom*Math.min(c.width,c.height)*.37/(1+depth*.18);return [c.width/2+x*f,c.height/2-z*f,depth];}
function line(ctx,points,c,stroke,width=1){if(points.length<2)return;ctx.strokeStyle=stroke;ctx.lineWidth=width;ctx.beginPath();points.forEach((p,i)=>{const s=project(p,c);if(i)ctx.lineTo(s[0],s[1]);else ctx.moveTo(s[0],s[1]);});ctx.stroke();}
function draw(){const f=DATA.frames[frame],bmax=shared?DATA.bmax:f.bmax,jmax=shared?DATA.jmax[component]:f.jmax[component];
 canvases.forEach((c,index)=>{const ratio=window.devicePixelRatio||1,w=Math.round(c.clientWidth*ratio),h=Math.round(c.clientHeight*ratio);if(c.width!==w||c.height!==h){c.width=w;c.height=h;}const ctx=c.getContext('2d');ctx.clearRect(0,0,w,h);ctx.fillStyle='#0e202d';ctx.fillRect(0,0,w,h);
  f.outline.forEach(p=>line(ctx,p,c,'#527080',.7*ratio));
  (DATA.chamber||[]).forEach(p=>line(ctx,p,c,'#7c8a94',.6*ratio));
  (DATA.coils||[]).forEach(p=>line(ctx,p,c,'#db9160',1.8*ratio));
  if(index===0){f.lines.forEach(path=>{for(let i=1;i<path.points.length;i++)line(ctx,[path.points[i-1],path.points[i]],c,color(bmax>0?path.field[i]/bmax:0),1.2*ratio);});}
  else{const points=f.current.map(p=>({p:project(p.slice(0,3),c),j:p[3+component]})).sort((a,b)=>b.p[2]-a.p[2]);ctx.globalAlpha=.52;for(const p of points){const scaled=jmax>0?p.j/jmax:0,t=component?.5+.5*scaled:scaled;ctx.fillStyle=color(t,component!==0);ctx.beginPath();ctx.arc(p.p[0],p.p[1],1.9*ratio,0,2*Math.PI);ctx.fill();}ctx.globalAlpha=1;}
 });
 ids('time').value=frame;ids('time-label').textContent=f.time_us.toFixed(3)+' µs';ids('metric-time').textContent=f.time_us.toFixed(3)+' µs';ids('field-max').textContent=(bmax*1e6).toPrecision(3)+' µT';
 const units=jmax>=1000?'kA/m²':'A/m²',scale=jmax>=1000?1000:1;ids('current-min').textContent=(component?-jmax/scale:0).toPrecision(3)+' '+units;ids('current-max').textContent=(jmax/scale).toPrecision(3)+' '+units;
 ids('current-title').textContent=['Current density magnitude','Radial current density','Azimuthal current density','Axial current density'][component];ids('current-gradient').style.background=component?'linear-gradient(90deg,#315a9e,#f1f5f5,#d35a47)':'linear-gradient(90deg,#3478ba,#67ddce,#fff1ad)';
 ids('metric-field').textContent=(f.bmax*1e6).toPrecision(4)+' µT';ids('metric-current').textContent=f.jmax[0].toPrecision(4)+' A/m²';ids('metric-temp').textContent=f.temperature.map(v=>v.toFixed(3)).join(' / ')+' eV';ids('scale').textContent=shared?'Use per-state color scale':'Use shared color scale';ids('status').textContent='Saved state '+(frame+1)+' of '+DATA.frames.length+' · '+(shared?'shared':'per-state')+' color scale · '+(DATA.status_note||'initialized flux, no formation claim');
 if(ids('drive'))ids('drive').textContent=f.drive_center_t===undefined?'':'Prescribed vacuum-centre field: '+(f.drive_center_t*1e3).toFixed(3)+' mT';}
ids('time').max=DATA.frames.length-1;ids('time').addEventListener('input',e=>{frame=+e.target.value;draw();});ids('component').addEventListener('change',e=>{component=+e.target.value;draw();});ids('scale').onclick=()=>{shared=!shared;draw();};ids('reset').onclick=()=>{yaw=-.65;pitch=.35;zoom=DATA.initial_zoom||1;draw();};
ids('play').onclick=()=>{if(playing){clearInterval(playing);playing=null;ids('play').textContent='Play saved states';}else{ids('play').textContent='Pause';playing=setInterval(()=>{frame=(frame+1)%DATA.frames.length;draw();},1000);}};
canvases.forEach(c=>{let drag=null;c.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY];c.setPointerCapture(e.pointerId);});c.addEventListener('pointermove',e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.008;pitch=Math.max(-1.45,Math.min(1.45,pitch+(e.clientY-drag[1])*.008));drag=[e.clientX,e.clientY];draw();});c.addEventListener('pointerup',()=>drag=null);c.addEventListener('pointercancel',()=>drag=null);c.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.4,Math.min(3,zoom*Math.exp(-e.deltaY*.001)));draw();},{passive:false});});
ids('scope').textContent=DATA.scope;window.addEventListener('resize',draw);draw();
