// Saved physical positions rendered through an orbiting perspective camera.
const CYAN='#67dccc', AMBER='#f2b56d';
export class Scene {
  constructor(canvas){
    this.canvas=canvas; this.ctx=canvas.getContext('2d'); this.reset(); this.count=1200;this.paths=[];this.samples=[];this.fieldMaximum=1;
    let previous=null;
    canvas.addEventListener('pointerdown',e=>{previous=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId)});
    canvas.addEventListener('pointermove',e=>{if(!previous)return;this.yaw+=(e.clientX-previous[0])*.008;this.pitch=Math.max(-1.4,Math.min(1.4,this.pitch+(e.clientY-previous[1])*.008));previous=[e.clientX,e.clientY];this.draw()});
    canvas.addEventListener('pointerup',()=>previous=null);canvas.addEventListener('pointercancel',()=>previous=null);
    canvas.addEventListener('wheel',e=>{e.preventDefault();this.zoom=Math.max(.5,Math.min(4,this.zoom*Math.exp(-e.deltaY*.001)));this.draw()},{passive:false});
    new ResizeObserver(()=>this.draw()).observe(canvas);
  }
  reset(){this.yaw=-.48;this.pitch=.22;this.zoom=1}
  project(p,w,h){
    const [x,y,z]=p,cy=Math.cos(this.yaw),sy=Math.sin(this.yaw),cp=Math.cos(this.pitch),sp=Math.sin(this.pitch);
    const xx=cy*x+sy*z, zz=-sy*x+cy*z, yy=cp*y-sp*zz, depth=sp*y+cp*zz;
    const scale=Math.min(w,h)*1.9*this.zoom/(1+depth/.8);
    return [w/2+xx*scale,h/2-yy*scale,depth,scale];
  }
  draw(){
    const [ctx,w,h]=canvasContext(this.canvas);ctx.clearRect(0,0,w,h);
    if(!this.geometry)return;
    const line=(points,color,width=1)=>{ctx.beginPath();points.forEach((p,i)=>{const q=this.project(p,w,h);i?ctx.lineTo(q[0],q[1]):ctx.moveTo(q[0],q[1])});ctx.strokeStyle=color;ctx.lineWidth=width;ctx.stroke()};
    const R=this.geometry.wall_radius;
    for(let k=0;k<3;k++){const pts=[];for(let i=0;i<=100;i++){const t=i*2*Math.PI/100,p=[R*Math.cos(t),R*Math.sin(t),0];if(k===1)pts.push([p[0],0,p[1]]);else if(k===2)pts.push([0,p[0],p[1]]);else pts.push(p)}line(pts,'#38536160')}
    const {vertices,edges}=this.geometry;
    for(const [a,b] of edges)line([vertices[a],vertices[b]],'#9bb2be88',1.6);
    for(const path of this.paths){
      for(let i=1;i<path.points.length;i+=3){
        const points=path.points.slice(Math.max(0,i-1),Math.min(path.points.length,i+3));
        const f=Math.log1p(9*path.values[i]/Math.max(this.fieldMaximum,1e-300))/Math.log(10);
        line(points,color(f),1.15);
      }
      const k=Math.floor(path.points.length*.6);
      if(k>0){const a=this.project(path.points[k-1],w,h),b=this.project(path.points[k],w,h),angle=Math.atan2(b[1]-a[1],b[0]-a[0]);ctx.strokeStyle='#c6e9e6aa';ctx.beginPath();ctx.moveTo(b[0]-5*Math.cos(angle-.5),b[1]-5*Math.sin(angle-.5));ctx.lineTo(b[0],b[1]);ctx.lineTo(b[0]-5*Math.cos(angle+.5),b[1]-5*Math.sin(angle+.5));ctx.stroke()}
    }
    const stride=Math.max(1,Math.ceil(this.samples.length/1200));
    for(let i=0;i<this.samples.length;i+=stride){const s=this.samples[i],q=this.project(s.point,w,h),norm=Math.hypot(...s.vector),f=Math.log1p(9*norm/Math.max(this.fieldMaximum,1e-300))/Math.log(10);ctx.fillStyle=color(f);ctx.globalAlpha=.75;ctx.beginPath();ctx.arc(q[0],q[1],2,0,2*Math.PI);ctx.fill();if(norm>0)line([s.point,s.point.map((x,k)=>x+.006*s.vector[k]/norm)],color(f),1)}ctx.globalAlpha=1;
    if(this.frame){
      // Stratified subset keeps both species visible without changing physical counts.
      const lists=[[],[]];this.frame.species.forEach((s,i)=>lists[s].push(i));
      const selected=[];const total=lists[0].length+lists[1].length;
      for(let s=0;s<2;s++){
        const n=Math.min(lists[s].length,Math.round(this.count*lists[s].length/Math.max(1,total)));
        for(let j=0;j<n;j++){const i=lists[s][Math.floor(j*lists[s].length/n)];selected.push({q:this.project(this.frame.positions_m[i],w,h),s})}
      }
      selected.sort((a,b)=>b.q[2]-a.q[2]);
      for(const {q,s} of selected){ctx.globalAlpha=s===0?.63:.85;ctx.fillStyle=s===0?CYAN:AMBER;ctx.beginPath();ctx.arc(q[0],q[1],s===0?1.15:1.65,0,2*Math.PI);ctx.fill()}
      ctx.globalAlpha=1;this.visible=selected.length;
    }
    const origin=this.project([-.15,-.13,0],w,h),axes=[[.024,0,0],[0,.024,0],[0,0,.024]];
    ctx.font='9px monospace';axes.forEach((p,i)=>{const end=this.project([p[0]-.15,p[1]-.13,p[2]],w,h);ctx.strokeStyle=['#bb7d77','#7fa994','#799dbb'][i];ctx.beginPath();ctx.moveTo(origin[0],origin[1]);ctx.lineTo(end[0],end[1]);ctx.stroke();ctx.fillStyle=ctx.strokeStyle;ctx.fillText(['x','y','z'][i],end[0]+3,end[1])});
  }
}

export function canvasContext(c){const w=c.clientWidth,h=c.clientHeight,dpr=Math.min(2,devicePixelRatio||1),pw=Math.round(w*dpr),ph=Math.round(h*dpr);if(c.width!==pw)c.width=pw;if(c.height!==ph)c.height=ph;const ctx=c.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);return [ctx,w,h]}
const stops=[[23,44,129],[25,125,172],[84,213,173],[237,230,139],[242,164,100]];
export function color(f){f=Math.max(0,Math.min(1,f))*(stops.length-1);const i=Math.min(stops.length-2,Math.floor(f)),t=f-i;return `rgb(${stops[i].map((v,k)=>Math.round(v*(1-t)+stops[i+1][k]*t)).join(',')})`}
export function drawSlice(canvas,frame,geometry){
  const [ctx,w,h]=canvasContext(canvas);ctx.clearRect(0,0,w,h);
  const values=frame.potential_slice_v,finite=values.filter(x=>x!==null),min=Math.min(0,...finite),max=Math.max(1e-12,...finite);
  const n=geometry.slice_axis_m.length,size=Math.min(w-70,h-24),left=(w-size)/2,top=8,cell=size/n;
  for(let iz=0;iz<n;iz++)for(let ix=0;ix<n;ix++){const v=values[iz*n+ix];if(v===null)continue;ctx.fillStyle=color((v-min)/(max-min));ctx.fillRect(left+ix*cell,top+(n-1-iz)*cell,cell+.3,cell+.3)}
  ctx.strokeStyle='#90aab55a';ctx.lineWidth=.5;ctx.beginPath();ctx.moveTo(left,top+size/2);ctx.lineTo(left+size,top+size/2);ctx.moveTo(left+size/2,top);ctx.lineTo(left+size/2,top+size);ctx.stroke();
  ctx.fillStyle='#879ca9';ctx.font='9px monospace';ctx.fillText('−19 cm',left,top+size+13);ctx.fillText('x',left+size/2,top+size+13);ctx.fillText('+19 cm',left+size-35,top+size+13);ctx.fillText('z',left-15,top+size/2);
  return [min,max];
}
export function drawHistory(canvas,history,index,target_ns=0){
  const [ctx,w,h]=canvasContext(canvas);ctx.clearRect(0,0,w,h);if(history.length<2)return;
  const keys=['electron_kinetic_j','ion_kinetic_j','field_j'],colors=[CYAN,AMBER,'#b1bfcc'];
  const maximum=Math.max(...history.flatMap(r=>keys.map(k=>r[k])),1e-30),end=history.at(-1).time_ns;
  const x=t=>42+(w-59)*t/end,y=v=>h-23-(h-43)*v/maximum;
  const unit=Math.max(end,target_ns)>=1000?'µs':'ns',scale=unit==='µs'?1000:1;
  ctx.font='8px monospace';ctx.fillStyle='#718a9a';ctx.fillText(`${(maximum*1e9).toPrecision(2)} nJ`,6,18);ctx.fillText('0',23,h-21);ctx.fillText(`0 ${unit}`,42,h-5);ctx.fillText(`${(end/scale).toPrecision(3)} ${unit}`,w-75,h-5);
  ctx.strokeStyle='#293b47';ctx.lineWidth=.5;for(let j=0;j<3;j++){const yy=20+j*(h-43)/2;ctx.beginPath();ctx.moveTo(42,yy);ctx.lineTo(w-17,yy);ctx.stroke()}
  keys.forEach((key,k)=>{ctx.strokeStyle=colors[k];ctx.lineWidth=1.4;ctx.beginPath();history.forEach((r,i)=>{i?ctx.lineTo(x(r.time_ns),y(r[key])):ctx.moveTo(x(r.time_ns),y(r[key]))});ctx.stroke()});
  ctx.strokeStyle='#e0ebed66';ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(x(history[index].time_ns),17);ctx.lineTo(x(history[index].time_ns),h-23);ctx.stroke();ctx.setLineDash([]);
}

// Actual saved radius and applied voltage; neither trace drives the simulation.
export function drawDrive(canvas,history,index,ceiling){
  const [ctx,w,h]=canvasContext(canvas);ctx.clearRect(0,0,w,h);
  if(history.length<2)return;
  const end=history.at(-1).time_ns/1000;
  const radii=history.map(r=>r.ion_motion?.radius_rms_m).filter(v=>v!==null&&Number.isFinite(v));
  const maxR=Math.max(...radii,1e-9)*1.08,maxV=Math.max(ceiling,1);
  const x=t=>58+(w-118)*t/end,y=f=>h-26-(h-50)*f;
  ctx.font='9px monospace';ctx.fillStyle='#f2b56d';ctx.fillText(`${(maxR*1000).toPrecision(3)} mm`,8,17);
  ctx.fillStyle='#67dccc';ctx.fillText(`${(maxV/1000).toPrecision(3)} kV`,w-59,17);
  ctx.fillStyle='#8395a5';ctx.fillText('0 µs',58,h-7);ctx.fillText(`${end.toPrecision(3)} µs`,w-105,h-7);
  ctx.strokeStyle='#293b47';ctx.lineWidth=.5;
  for(let j=0;j<3;j++){ctx.beginPath();ctx.moveTo(58,y(j/2));ctx.lineTo(w-60,y(j/2));ctx.stroke()}
  for(const [color,value] of [['#f2b56d',r=>r.ion_motion?.radius_rms_m==null?NaN:r.ion_motion.radius_rms_m/maxR],['#67dccc',r=>r.voltage_v/maxV]]){
    ctx.strokeStyle=color;ctx.lineWidth=1.5;ctx.beginPath();let connected=false;
    for(const r of history){const v=value(r);if(!Number.isFinite(v)){connected=false;continue}const xx=x(r.time_ns/1000),yy=y(v);connected?ctx.lineTo(xx,yy):ctx.moveTo(xx,yy);connected=true}ctx.stroke();
  }
  ctx.strokeStyle='#e0ebed66';ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(x(history[index].time_ns/1000),20);ctx.lineTo(x(history[index].time_ns/1000),h-26);ctx.stroke();
}
