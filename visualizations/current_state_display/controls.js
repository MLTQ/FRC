// Camera presets and keyboard access; no changes to archived field data.
const presets={perspective:{yaw:.50,pitch:.30,distance:3.70},side:{yaw:0,pitch:0,distance:3.70},end:{yaw:Math.PI/2,pitch:0,distance:3.70}};
for(const button of document.querySelectorAll('[data-view]'))button.onclick=()=>{
  Object.assign(state,presets[button.dataset.view],{orbit:false});
  document.getElementById('orbit').setAttribute('aria-pressed','false');
  dirty=true;
};
for(const canvas of document.querySelectorAll('canvas')){
  canvas.tabIndex=0;
  canvas.setAttribute('aria-description','Use arrow keys to rotate and plus or minus to zoom. Both panels move together.');
  canvas.addEventListener('keydown',e=>{
    if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','=','-'].includes(e.key))return;
    e.preventDefault();state.orbit=false;
    document.getElementById('orbit').setAttribute('aria-pressed','false');
    if(e.key==='ArrowLeft')state.yaw-=.10;
    if(e.key==='ArrowRight')state.yaw+=.10;
    if(e.key==='ArrowUp')state.pitch=Math.min(1.25,state.pitch+.10);
    if(e.key==='ArrowDown')state.pitch=Math.max(-1.25,state.pitch-.10);
    if(e.key==='+'||e.key==='=')state.distance=Math.max(2.8,state.distance/1.1);
    if(e.key==='-')state.distance=Math.min(7,state.distance*1.1);
    dirty=true;
  });
}
