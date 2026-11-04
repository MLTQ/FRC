import {Scene,drawSlice,drawHistory,drawDrive} from './scene.js?v=20261027-4';
import {FieldsView} from './fields.js?v=20261027-4';
const $=s=>document.querySelector(s),scene=new Scene($('#scene'));
const fields=new FieldsView();
let state=null,experiment='iec',selected=null,geometry=null,history=[],index=0,frame=null,playing=false,selection=0,frameRequest=0,loading=false,frameLoading=false;
const frameCache=new Map();
fields.onModeChange=()=>{if(selected&&history.length)showFrame(+$('#time').value)};
const number=(v,d=2)=>Number(v).toLocaleString(undefined,{maximumFractionDigits:d});
const timeLabel=(ns,target=ns)=>target>=1000?`${number(ns/1000,3)} µs`:`${number(ns,2)} ns`;
const runLabel=r=>(r.label||'').replace('HOPS · hops','HOPS · radius feedback');
const statusLabel=s=>({completed:'Completed',running:'Running',initializing:'Preparing field',meshing:'Meshing',failed:'Failed',cancelled:'Stopped',saved:'Saved'})[s]||s;
function note(text){$('#notice').hidden=!text;$('#notice').textContent=text||''}
async function get(url){const response=await fetch(url,{cache:'no-store'});if(!response.ok)throw Error(`Could not read saved data (${response.status}).`);return response.json()}
async function post(url,values){const response=await fetch(url,{method:'POST',headers:{'Content-Type':'application/json','X-Lab-Token':state.token},body:JSON.stringify(values)});const data=await response.json();if(!response.ok)throw Error(data.error||'Request failed');return data}
function renderNavigation(){
  $('#experiments').replaceChildren(...state.experiments.map(e=>{const b=document.createElement('button');b.className='experiment'+(e.id===experiment?' active':'');b.innerHTML=`<span class="icon">${e.id==='iec'?'◈':'◎'}</span><span><strong></strong><small></small></span>`;b.querySelector('strong').textContent=e.name;b.querySelector('small').textContent=e.family;b.onclick=()=>chooseExperiment(e.id);return b}));
  const rows=state.runs.filter(r=>r.experiment===experiment);$('#run-count').textContent=rows.length;
  $('#runs').replaceChildren(...rows.map(r=>{const b=document.createElement('button');b.className='run'+(selected?.id===r.id?' active':'');const strong=document.createElement('strong');strong.textContent=r.config?.drive_mode&&r.config.drive_mode!=='legacy'?r.label:r.config?`${number(r.config.bias_v)} V · ${r.config.source_ua>0?'source '+number(r.config.source_ua)+' µA':'initial plasma'}`:r.label;const small=document.createElement('small');const dot=document.createElement('span');dot.className='dot '+r.status;const detail=r.config?`${timeLabel(r.target_ns)} · ${number(r.config.dt_ps)} ps`:`${number(r.frames||0,0)} states`;small.append(dot,document.createTextNode(`${statusLabel(r.status)} · ${detail}`));b.append(strong,small);b.onclick=()=>chooseRun(r);return b}));
}
async function chooseExperiment(id,preferred=null){
  experiment=id;selected=null;selection++;playing=false;$('#play').textContent='▶';
  const e=state.experiments.find(e=>e.id===id);$('#family').textContent=id==='iec'?'ELECTROSTATIC':'MAGNETIC';$('#title').textContent=e.name;$('#description').textContent=e.description;$('#model-tag').textContent=id==='iec'?'3D KINETIC MODEL':'FORMATION & PARTICLE MODELS';$('#configure').hidden=id!=='iec';$('#scope').textContent=e.scope;
  const rows=state.runs.filter(r=>r.experiment===id);renderNavigation();$('#empty').hidden=rows.length>0||id!=='iec';$('#iec-workspace').hidden=true;$('#frc-workspace').hidden=id!=='frc';
  if(rows.length)await chooseRun(rows.find(r=>r.id===preferred)||rows[0]);
}
async function chooseRun(run){
  const ticket=++selection;selected=run;fields.clear();frameCache.clear();frame=null;geometry=null;history=[];scene.geometry=null;scene.frame=null;scene.draw();index=0;playing=false;$('#play').textContent='▶';$('#metrics').replaceChildren();note('');renderNavigation();
  $('#time').max=0;$('#time').value=0;$('#time-label').textContent='0 ns';$('#frame-label').textContent='Preparing saved states…';$('#count-label').textContent='0';$('#balance').textContent='';
  for(const id of ['slice','history','drive-history']){const c=$('#'+id);c.getContext('2d').clearRect(0,0,c.width,c.height)}
  $('#hops-panel').hidden=true;
  $('#slice-min').textContent='';$('#slice-max').textContent='';
  $('#empty').hidden=true;$('#iec-workspace').hidden=run.experiment!=='iec';$('#frc-workspace').hidden=run.experiment!=='frc';
  window.history.replaceState(null,'','?run='+encodeURIComponent(run.id));
  if(run.experiment==='frc'){$('#frc-view').src=run.viewer_url;return}
  renderRun();
  try{
    if(run.frames){const data=await Promise.all([get(run.url+'/geometry.json'),get(run.url+'/history.json')]);if(ticket!==selection)return;[geometry,history]=data;scene.geometry=geometry;$('#time').max=history.length-1;index=history.length-1;await showFrame(index,ticket)}
  }catch(error){if(ticket===selection)note(error.message)}
}
function renderRun(){
  if(!selected||selected.experiment!=='iec')return;
  $('#run-state').textContent=statusLabel(selected.status);$('#run-name').textContent=runLabel(selected);$('#run-id').textContent=selected.id;
  const isHops=selected.config?.drive_mode&&selected.config.drive_mode!=='legacy';$('#title').textContent=isHops?'HOPS':'Positive cage';$('#description').textContent=isHops?'Harmonically Oscillating Plasma Sphere · positive cage, adaptive timing and matched controls.':state.experiments.find(e=>e.id==='iec').description;
  $('#stop').hidden=!['running','initializing','meshing'].includes(selected.status)||(selected.config?.drive_mode&&selected.config.drive_mode!=='legacy');$('#end-label').textContent=`${timeLabel(selected.target_ns)} target`;
  const r=selected.resolution;$('#resolution').textContent=geometry?`${number(geometry.nodes,0)} nodes · ${number(geometry.tetrahedra,0)} tetrahedra · C = ${number(geometry.capacitance_pf,2)} pF`:'Preparing the conforming 3D mesh…';
  if(r)$('#resolution').textContent+=`\nInitial Debye length ${number(r.initial_debye_m*1000,1)} mm · ωpe Δt = ${number(r.initial_omega_pe_dt,4)}`;
  const d=frame?.diagnostics,particleScale=d?d.electron_kinetic_j+d.ion_kinetic_j+d.absorbed_j+d.injected_j:0;
  if(particleScale>0)$('#resolution').textContent+=`\nBalance error / particle-energy throughput: ${number(100*Math.abs(d.energy_residual_j)/particleScale,4)}%`;
  const config=selected.config||{},controlled=config.drive_mode&&config.drive_mode!=='legacy';
  const parameters=state.experiments.find(e=>e.id==='iec').parameters.filter(p=>!controlled||!['rise_ns','modulation','frequency_mhz'].includes(p.key));
  if(controlled)parameters.push(...[['warmup_ns','Common warmup','ns'],['response_ns','Voltage response','ns'],['sensor_ns','Radius sensor','ns'],['release_fraction','Release radius / reference',''],['return_fraction','Return radius / reference','']].map(([key,label,unit])=>({key,label,unit})));
  $('#saved-inputs').replaceChildren(...parameters.flatMap(p=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=p.label;const v=config[p.key];dd.textContent=`${Math.abs(v)>=1e7?Number(v).toExponential(2):number(v,4)} ${p.unit}`;return [dt,dd]}));
  if(selected.error)note(selected.error);
}
async function showFrame(i,ticket=selection){
  if(!selected||!history.length)return;
  const request=++frameRequest;
  frameLoading=true;
  i=Math.max(0,Math.min(history.length-1,Math.round(i)));$('#time').value=i;
  const url=selected.url+`/frame_${String(i).padStart(4,'0')}.json`;
  try{
    const next=frameCache.get(url)||await get(url);
    if(ticket!==selection||request!==frameRequest)return;
    const prepared=await fields.prepare(selected,i,geometry,history);
    if(!prepared||ticket!==selection||request!==frameRequest)return;
    frameCache.set(url,next);while(frameCache.size>8)frameCache.delete(frameCache.keys().next().value);
    await new Promise(resolve=>requestAnimationFrame(()=>{
      if(ticket===selection&&request===frameRequest&&fields.commit(prepared)){index=i;frame=next;scene.frame=frame;draw()}
      resolve();
    }));
  }catch(error){if(ticket===selection)note(error.message)}finally{if(request===frameRequest)frameLoading=false}
}
function draw(){
  if(!frame||!geometry)return;scene.count=+$('#count').value;scene.draw();$('#count-label').textContent=number(scene.visible,0);
  const d=frame.diagnostics;$('#time-label').textContent=timeLabel(d.time_ns,selected.target_ns);$('#frame-label').textContent=`State ${index+1} / ${history.length}`;
  const initialIons=history[0]?.ions||d.ions;
  const data=[['Ions retained',initialIons?100*d.ions/initialIons:0,'%'],['Core ions · r < 2 cm',d.core_ions,'macro'],['Cage bias',d.voltage_v,'V'],['Centre potential',d.centre_potential_v,'V']];
  if(d.ion_motion?.radius_rms_m!==null&&d.ion_motion){data.push(['Ion RMS radius',d.ion_motion.radius_rms_m*100,'cm'],['Radial speed',d.ion_motion.radial_speed_m_s,'m/s'],['Breathing energy',d.ion_motion.breathing_j*1e9,'nJ'],['Residual ion motion',d.ion_motion.residual_motion_j*1e9,'nJ'])}
  $('#metrics').replaceChildren(...data.map(([label,value,unit])=>{const div=document.createElement('div');div.className='metric';const small=document.createElement('small');small.textContent=label;const strong=document.createElement('strong');strong.textContent=value!==0&&Math.abs(value)<.01?value.toExponential(2):number(value,unit==='macro'?0:2);const span=document.createElement('span');span.textContent=unit;strong.append(span);div.append(small,strong);return div}));
  const [min,max]=drawSlice($('#slice'),frame,geometry);$('#slice-min').textContent=`${number(min,1)} V`;$('#slice-max').textContent=`${number(max,1)} V`;
  const controlled=selected.config.drive_mode&&selected.config.drive_mode!=='legacy';$('#hops-panel').hidden=!controlled;if(controlled)drawDrive($('#drive-history'),history,index,selected.config.bias_v);
  drawHistory($('#history'),history,index,selected.target_ns);const residual=100*d.energy_residual_fraction;$('#balance').textContent=`Residual ${residual!==0&&Math.abs(residual)<.001?residual.toExponential(1):number(residual,3)}%`;
  renderRun();
}
async function refresh(){
  if(loading)return;loading=true;
  try{state=await get('/api/state');renderNavigation();if(!selected){const requested=state.runs.find(r=>r.id===new URLSearchParams(location.search).get('run'));await chooseExperiment(requested?.experiment||experiment,requested?.id);return}
    const current=state.runs.find(r=>r.id===selected.id);if(current&&current.experiment==='iec'){
      const old=selected.frames;selected=current;renderRun();if(current.frames>old){const atEnd=index>=history.length-1;if(!geometry){geometry=await get(current.url+'/geometry.json');scene.geometry=geometry}history=await get(current.url+'/history.json');$('#time').max=history.length-1;if(atEnd)await showFrame(history.length-1)}
    }
  }catch(error){note(error.message)}finally{loading=false}
}
function setup(){
  const e=state.experiments.find(e=>e.id==='iec');$('#parameters').replaceChildren(...e.parameters.map(p=>{const label=document.createElement('label');label.textContent=p.label;const unit=document.createElement('span');unit.textContent=p.unit;const input=document.createElement('input');input.type='number';input.name=p.key;input.value=p.default;input.min=p.min;input.max=p.max;input.step='any';input.required=true;label.append(unit,input);return label}));$('#form-error').textContent='';$('#setup').showModal();
}
$('#configure').onclick=setup;$('#empty-start').onclick=setup;$('#close-setup').onclick=()=>$('#setup').close();$('#cancel-setup').onclick=()=>$('#setup').close();$('#refresh').onclick=refresh;
$('#setup-form').onsubmit=async e=>{e.preventDefault();$('#launch').disabled=true;$('#form-error').textContent='';try{const values=Object.fromEntries([...new FormData(e.target)].map(([k,v])=>[k,Number(v)]));const result=await post('/api/launch',values);$('#setup').close();note('Simulation started. Preparing the mesh and initial field.');for(let attempt=0;attempt<10;attempt++){await new Promise(r=>setTimeout(r,500));state=await get('/api/state');const run=state.runs.find(r=>r.id===result.id);if(run){await chooseRun(run);break}}}catch(error){$('#form-error').textContent=error.message}finally{$('#launch').disabled=false}};
$('#stop').onclick=async()=>{try{await post('/api/stop',{id:selected.id});note('Stop requested. The run will retain its saved states.')}catch(error){note(error.message)}};
$('#time').oninput=()=>{playing=false;$('#play').textContent='▶';showFrame(+$('#time').value)};$('#count').oninput=draw;$('#reset-camera').onclick=()=>{scene.reset();scene.draw()};
$('#play').onclick=()=>{playing=!playing;$('#play').textContent=playing?'Ⅱ':'▶';if(playing){const target=+$('#time').value;if(target>=history.length-1)showFrame(0);else if(target!==index)showFrame(target)}};
$('#frc-view').onload=()=>{const control=$('#frc-view').contentDocument?.querySelector('#run');control?.addEventListener('change',()=>{const run=state.runs.find(r=>r.experiment==='frc'&&r.id===control.value);if(run){selected=run;renderNavigation();window.history.replaceState(null,'','?run='+encodeURIComponent(run.id))}})};
setInterval(()=>{if(playing&&history.length&&!fields.busy&&!frameLoading){if(index>=history.length-1){playing=false;$('#play').textContent='▶'}else showFrame(index+1)}},140);
window.addEventListener('resize',()=>{if(frame)draw()});
await refresh();setInterval(refresh,5000);
