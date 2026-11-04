import {Scene,canvasContext,color} from './scene.js?v=20261027-4';
const $=s=>document.querySelector(s),cyan='#67dccc',amber='#f2b56d';
const fmt=x=>!Number.isFinite(x)?'—':x===0?'0':Math.abs(x)<.01||Math.abs(x)>=1e5?x.toExponential(2):x.toLocaleString(undefined,{maximumFractionDigits:2});
const timeUs=r=>r.time_us??r.time_ns/1000;
const timeText=t=>Number(t).toLocaleString(undefined,{maximumFractionDigits:6});
const quantities={j:['Total |J|','A/m²'],je:['Electron |Jₑ|','A/m²'],ji:['Ion |Jᵢ|','A/m²'],ne:['Electron density','m⁻³'],ni:['Ion density','m⁻³'],rho:['Charge density','C/m³'],electric:['Electric |E|','V/m']};

function axes(canvas,xmax,ymin,ymax,xlabel,ylabel){
  const [ctx,w,h]=canvasContext(canvas);ctx.clearRect(0,0,w,h);
  const x=t=>55+(w-75)*t/Math.max(xmax,1e-30),y=t=>h-36-(h-65)*(t-ymin)/Math.max(ymax-ymin,1e-30);
  ctx.font='10px monospace';ctx.fillStyle='#8ea5b3';ctx.fillText(ylabel,12,17);ctx.fillText(xlabel,w/2-30,h-7);
  ctx.strokeStyle='#2b3f4c';ctx.lineWidth=.6;
  for(let i=0;i<4;i++){const value=ymin+(ymax-ymin)*i/3,yy=y(value);ctx.beginPath();ctx.moveTo(55,yy);ctx.lineTo(w-20,yy);ctx.stroke();ctx.fillText(fmt(value),4,yy+3)}
  for(let i=0;i<4;i++){const value=xmax*i/3;ctx.fillText(fmt(value),x(value)-10,h-21)}
  return {ctx,w,h,x,y};
}

function drawPhase(canvas,points){
  const maxR=Math.max(.05,...points.map(p=>p[0]))*100,maxV=Math.max(.01,...points.map(p=>Math.abs(p[1])/1000));
  const {ctx,x,y}=axes(canvas,maxR,-maxV*1.1,maxV*1.1,'Radius [cm]','Radial velocity [km/s]');
  ctx.strokeStyle='#d9e3e955';ctx.beginPath();ctx.moveTo(x(0),y(0));ctx.lineTo(x(maxR),y(0));ctx.stroke();
  ctx.fillStyle=amber;ctx.globalAlpha=.65;
  for(const [r,v] of points){ctx.beginPath();ctx.arc(x(r*100),y(v/1000),1.8,0,Math.PI*2);ctx.fill()}ctx.globalAlpha=1;
}

function drawProfile(canvas,profile){
  const maximum=Math.max(1,...profile.ne,...profile.ni),logMax=Math.log10(1+maximum);
  const {ctx,x,y}=axes(canvas,20,0,logMax,'Radius [cm]','log₁₀(1 + density [m⁻³])');
  for(const [key,stroke] of [['ne',cyan],['ni',amber]]){ctx.strokeStyle=stroke;ctx.lineWidth=1.6;ctx.beginPath();profile[key].forEach((v,i)=>{const xx=x(profile.radius_m[i]*100),yy=y(Math.log10(1+v));i?ctx.lineTo(xx,yy):ctx.moveTo(xx,yy)});ctx.stroke()}
}

function drawPopulation(canvas,history,index){
  const maximum=Math.max(1,...history.flatMap(r=>[r.electrons,r.ions])),end=timeUs(history.at(-1))||1;
  const {ctx,x,y,h}=axes(canvas,end,0,Math.log10(maximum+1),'Time [µs]','log₁₀(1 + live macros)');
  for(const [key,stroke] of [['electrons',cyan],['ions',amber]]){ctx.strokeStyle=stroke;ctx.lineWidth=1.7;ctx.beginPath();history.forEach((r,i)=>{i?ctx.lineTo(x(timeUs(r)),y(Math.log10(1+r[key]))):ctx.moveTo(x(timeUs(r)),y(Math.log10(1+r[key])))});ctx.stroke()}
  ctx.strokeStyle='#d9e3e988';ctx.setLineDash([4,4]);ctx.beginPath();ctx.moveTo(x(timeUs(history[index])),25);ctx.lineTo(x(timeUs(history[index])),h-36);ctx.stroke();ctx.setLineDash([]);
}

function drawMoment(canvas,slice,key){
  const [ctx,w,h]=canvasContext(canvas);ctx.clearRect(0,0,w,h);
  const values=slice.values[key],finite=values.filter(v=>v!==null),max=Math.max(1e-300,...finite.map(Math.abs)),signed=key==='rho';
  const n=slice.axis_m.length,size=Math.min(w-60,h-50),left=(w-size)/2,top=10,cell=size/n;
  for(let iz=0;iz<n;iz++)for(let ix=0;ix<n;ix++){const v=values[iz*n+ix];if(v===null)continue;const f=signed?.5+.5*v/max:Math.log1p(99*Math.max(0,v)/max)/Math.log(100);ctx.fillStyle=color(f);ctx.fillRect(left+ix*cell,top+(n-1-iz)*cell,cell+.2,cell+.2)}
  ctx.font='10px monospace';ctx.fillStyle='#8ea5b3';ctx.fillText('−18 cm',left,top+size+17);ctx.fillText('x',left+size/2,top+size+17);ctx.fillText('+18 cm',left+size-40,top+size+17);ctx.fillText('z',left-15,top+size/2);
  return `${signed?'±':''}${fmt(max)} ${quantities[key][1]} · ${signed?'linear signed':'log color'} · sampled range`;
}

export class FieldsView {
  constructor(){
    this.scene=new Scene($('#field-scene'));this.scene.count=0;this.busy=false;this.cache=new Map();this.generation=0;this.request=0;
    $('#field-kind').onchange=()=>this.onModeChange?.();
    $('#moment-kind').onchange=()=>this.render();$('#phase-kind').onchange=()=>this.render();
    $('#field-reset').onclick=()=>{this.scene.reset();this.scene.draw()};
  }
  clear(){
    this.generation++;this.request++;this.params=null;this.data=null;
    if(this.wanted)this.wanted.resolve(null);this.wanted=null;this.shown=null;
    this.scene.geometry=null;this.scene.paths=[];this.scene.samples=[];this.scene.draw();
    for(const id of ['moment-slice','phase-space','population-history']){const c=$('#'+id);c.getContext('2d').clearRect(0,0,c.width,c.height)}
    $('#field-status').textContent='Choose a saved state to inspect its native fields.';$('#moment-scale').textContent='';
  }
  prepare(run,index,geometry,history){
    const mode=$('#field-kind').value,magnetic=mode==='magnetic',key=`${run.id}:${index}:${magnetic}`;
    const request=++this.request;
    if(this.wanted)this.wanted.resolve(null);
    return new Promise(resolve=>{
      this.wanted={run,index,geometry,history,mode,magnetic,key,generation:this.generation,request,resolve};
      if(!this.busy)this.pump();
    });
  }
  commit(item){
    if(item.generation!==this.generation||item.request!==this.request||item.mode!==$('#field-kind').value)return false;
    this.data=item.data;this.params=[item.run,item.index,item.geometry,item.history];this.shown=item.key;
    this.mode=item.mode;
    this.scene.geometry=item.geometry;this.render();return true;
  }
  async pump(){
    this.busy=true;
    try{while(this.wanted){
      const item=this.wanted;this.wanted=null;
      try{
        let data=this.cache.get(item.key);
        if(!data){
          if(!this.data)$('#field-status').textContent=`Reading full saved state at ${timeText(timeUs(item.history[item.index]))} µs…`;
          const r=await fetch(`/api/observables?run=${encodeURIComponent(item.run.id)}&frame=${item.index}&magnetic=${item.magnetic?1:0}`);data=await r.json();if(!r.ok)throw Error(data.error||'Native field data unavailable');this.cache.set(item.key,data);while(this.cache.size>8)this.cache.delete(this.cache.keys().next().value);
        }
        item.resolve(item.generation===this.generation&&item.request===this.request?{...item,data}:null);
      }catch(error){
        if(item.generation===this.generation&&item.request===this.request)$('#field-status').textContent=`${this.data?`Showing ${timeText(this.data.time_us)} µs · `:''}${error.message}`;
        item.resolve(null);
      }
    }}finally{this.busy=false;if(this.wanted)this.pump()}
  }
  render(){
    if(!this.data||!this.params)return;
    const d=this.data,mode=this.mode;
    this.scene.paths=[];this.scene.samples=[];
    if(mode==='current'){
      this.scene.samples=d.current.points.map((point,i)=>({point,vector:d.current.vectors[i]}));this.scene.fieldMaximum=d.current.max_am2;
      $('#field-status').textContent=`${timeText(d.time_us)} µs · total current moment · max ${fmt(d.current.max_am2)} A/m²`;
      $('#field-scale-max').textContent=`${fmt(d.current.max_am2)} A/m²`;
      $('#field-note').textContent='Current from all saved particles. Glyph length is normalized for display; color carries magnitude. No trajectory-current continuity claim.';
    }else if(mode==='magnetic'){
      if(!d.magnetic)return;
      this.scene.paths=d.magnetic.lines;this.scene.fieldMaximum=d.magnetic.max_t;
      $('#field-status').textContent=`${timeText(d.time_us)} µs · plasma-only B · max sampled ${fmt(d.magnetic.max_t*1e9)} nT`;
      $('#field-scale-max').textContent=`${fmt(d.magnetic.max_t*1e9)} nT`;
      $('#field-note').textContent=d.magnetic.scope+` Softening ${(d.magnetic.softening_m*1000).toFixed(0)} mm; ${d.magnetic.grid_points}³ display samples.`;
    }else{
      this.scene.paths=d.electric.lines;this.scene.fieldMaximum=Math.max(1e-300,...d.electric.lines.flatMap(p=>p.values));
      $('#field-status').textContent=`${timeText(d.time_us)} µs · E = −∇φ · native peak ${fmt(d.electric.max_vm)} V/m`;
      $('#field-scale-max').textContent=`${fmt(this.scene.fieldMaximum)} V/m · line samples`;
      $('#field-note').textContent='Lines follow the saved 3D FEM electric field and end at steel boundaries. Arrows point along E. Color shows |E|; neither seed density nor line count measures charge density.';
    }
    this.scene.draw();const key=$('#moment-kind').value;$('#moment-scale').textContent=drawMoment($('#moment-slice'),d.slice,key);
    if($('#phase-kind').value==='radial'){drawProfile($('#phase-space'),d.radial);$('#phase-note').textContent=d.radial.scope+' Cyan: electrons; amber: ions.'}
    else{drawPhase($('#phase-space'),d.phase);$('#phase-note').textContent=`All ${d.phase.length} live ions, measured about the ion centroid. Negative velocity is inward. Losses can bias the surviving cloud; these points are not particle tracks.`}
    drawPopulation($('#population-history'),this.params[3],this.params[1]);
  }
}
