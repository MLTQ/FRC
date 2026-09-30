/* Offline WebGL2 renderer for archived axisymmetric SI fields. */
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const unit=a=>a.map(v=>v/Math.hypot(...a));
const VS=`#version 300 es
in vec2 position;out vec2 uv;void main(){uv=position;gl_Position=vec4(position,0.,1.);}`;
const FS=`#version 300 es
precision highp float;in vec2 uv;out vec4 color;
uniform sampler2D field;uniform vec3 eye,right,up,forward;
uniform float aspect,tanHalf,jmax,density;uniform vec2 extent;uniform int cutaway;
vec4 lookup(vec2 p){vec2 s=p*vec2(textureSize(field,0))-.5;ivec2 i=ivec2(floor(s));vec2 f=fract(s);ivec2 hi=textureSize(field,0)-1;
return mix(mix(texelFetch(field,clamp(i,ivec2(0),hi),0),texelFetch(field,clamp(i+ivec2(1,0),ivec2(0),hi),0),f.x),
mix(texelFetch(field,clamp(i+ivec2(0,1),ivec2(0),hi),0),texelFetch(field,clamp(i+ivec2(1,1),ivec2(0),hi),0),f.x),f.y);}
vec3 heat(float q){vec3 a=vec3(.28,.015,.008),b=vec3(1.,.18,.025),c=vec3(1.,.72,.24),d=vec3(1.,.96,.77);
return q<.4?mix(a,b,q/.4):q<.78?mix(b,c,(q-.4)/.38):mix(c,d,(q-.78)/.22);}
void main(){vec3 rd=normalize(forward+tanHalf*(uv.x*aspect*right+uv.y*up));
vec3 bg=mix(vec3(.017,.028,.048),vec3(.032,.055,.079),exp(-dot(uv,uv)*1.4));
float b=dot(eye,rd),c=dot(eye,eye)-1.09*1.09,h=b*b-c;vec3 sum=vec3(0);float alpha=0.;
if(h>0.&&density>0.){float start=max(0.,-b-sqrt(h)),stop=-b+sqrt(h),ds=(stop-start)/240.;
for(int i=0;i<240;i++){vec3 p=eye+rd*(start+(float(i)+.5)*ds);if(cutaway==1&&p.z>0.&&p.y>0.)continue;
vec2 t=vec2(length(p.yz)/extent.x,(p.x/extent.y+1.)*.5);float q=clamp(lookup(t).r/jmax,0.,1.);
float a=1.-exp(-density*ds*6.*pow(q,1.65));sum+=(1.-alpha)*a*heat(q)*1.35;alpha+=(1.-alpha)*a;if(alpha>.995)break;}}
color=vec4(sum+bg*(1.-alpha),1.);}`;
const LINEVS=`#version 300 es
precision highp float;in vec3 position,normal;in float value;uniform mat4 camera;uniform float width;
out vec3 n;out vec3 world;out float v;void main(){world=position;n=normal;v=value;gl_Position=camera*vec4(position+normal*width,1.);}`;
const LINEFS=`#version 300 es
precision highp float;in vec3 n;in vec3 world;in float v;out vec4 color;uniform float opacity;uniform int mode,cutaway;
void main(){if(mode==0&&cutaway==1&&world.z>0.&&world.y>0.)discard;float q=clamp(v,0.,1.);vec3 c;
if(mode==0)c=mix(vec3(.97,.18,.02),vec3(1.,.93,.62),q);
else c=q<.4?mix(vec3(.06,.25,.65),vec3(.06,.91,.89),q/.4):mix(vec3(.06,.91,.89),vec3(1.,.91,.59),(q-.4)/.6);
float light=.55+.45*abs(dot(normalize(n),normalize(vec3(.3,.8,1.))));color=vec4(c*light,opacity);}`;
function program(gl,vs,fs){const compile=(t,s)=>{const x=gl.createShader(t);gl.shaderSource(x,s);gl.compileShader(x);if(!gl.getShaderParameter(x,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(x));return x;};
const p=gl.createProgram();gl.attachShader(p,compile(gl.VERTEX_SHADER,vs));gl.attachShader(p,compile(gl.FRAGMENT_SHADER,fs));gl.linkProgram(p);if(!gl.getProgramParameter(p,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(p));return p;}
function multiply(a,b){const c=new Float32Array(16);for(let j=0;j<4;j++)for(let i=0;i<4;i++)for(let k=0;k<4;k++)c[j*4+i]+=a[k*4+i]*b[j*4+k];return c;}
function camera(yaw,pitch,distance,aspect){const eye=[distance*Math.sin(yaw)*Math.cos(pitch),distance*Math.sin(pitch),distance*Math.cos(yaw)*Math.cos(pitch)];
const f=unit(eye.map(v=>-v)),r=unit(cross(f,[0,1,0])),u=cross(r,f);
const view=new Float32Array([r[0],u[0],-f[0],0,r[1],u[1],-f[1],0,r[2],u[2],-f[2],0,-dot(r,eye),-dot(u,eye),dot(f,eye),1]);
const t=Math.tan(36*Math.PI/360),near=.05,far=30;const perspective=new Float32Array([1/(aspect*t),0,0,0,0,1/t,0,0,0,0,(far+near)/(near-far),-1,0,0,2*far*near/(near-far),0]);
return {eye,r,u,f,t,matrix:multiply(perspective,view)};}
function tube(vertices,points,values){const sides=6;for(let i=1;i<points.length;i++){
const a=points[i-1],b=points[i],delta=b.map((x,k)=>x-a[k]);if(Math.hypot(...delta)<1e-9)continue;
const t=unit(delta),n=unit(cross(t,Math.abs(t[1])<.9?[0,1,0]:[1,0,0])),v=cross(t,n);
const normal=k=>n.map((x,j)=>x*Math.cos(k*2*Math.PI/sides)+v[j]*Math.sin(k*2*Math.PI/sides));
const put=(p,n,q)=>vertices.push(...p,...n,q);
for(let k=0;k<sides;k++){let n0=normal(k),n1=normal(k+1);put(a,n0,values[i-1]);put(b,n0,values[i]);put(b,n1,values[i]);put(a,n0,values[i-1]);put(b,n1,values[i]);put(a,n1,values[i-1]);}}}
class FieldView {
 constructor(canvas,data,mode){this.canvas=canvas;this.data=data;this.mode=mode;const gl=this.gl=canvas.getContext('webgl2',{antialias:true,preserveDrawingBuffer:true,alpha:false});if(!gl)throw Error('WebGL 2 is needed for this viewer.');
this.volume=program(gl,VS,FS);this.line=program(gl,LINEVS,LINEFS);this.quad=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,this.quad);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),gl.STATIC_DRAW);
this.buffer=gl.createBuffer();this.texture=gl.createTexture();this.setFrame(2);}
 setFrame(index){const {gl,data}=this;const f=data.frames[index];gl.bindTexture(gl.TEXTURE_2D,this.texture);
const bytes=Uint8Array.from(atob(f.texture),c=>c.charCodeAt(0));gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA32F,...data.texture_size,0,gl.RGBA,gl.FLOAT,new Float32Array(bytes.buffer));
gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);
const vertices=[];
if(this.mode===0){for(const [r,z,q] of f.rings){let p=[],v=[];for(let k=0;k<=100;k++){let phi=k*2*Math.PI/100;p.push([z,r*Math.cos(phi),r*Math.sin(phi)]);v.push(q);}tube(vertices,p,v);}}
else {for(const path of f.lines){let count=path.closed?5:3;for(let k=0;k<count;k++){let phi=k*2*Math.PI/count+.16;const points=path.points.map(([z,r])=>[z,r*Math.cos(phi),r*Math.sin(phi)]);tube(vertices,points,path.b.map(b=>b/data.bmax_mT));}}}
this.count=vertices.length/7;gl.bindBuffer(gl.ARRAY_BUFFER,this.buffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(vertices),gl.STATIC_DRAW);}
 draw(state,exportSize){const {gl,canvas}=this;const rect=canvas.getBoundingClientRect();const scale=Math.min(devicePixelRatio||1,1.6);
const w=exportSize?.[0]||Math.round(rect.width*scale),h=exportSize?.[1]||Math.round(rect.height*scale);if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}
gl.viewport(0,0,w,h);gl.clearColor(.02,.03,.05,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);const c=camera(state.yaw,state.pitch,state.distance,w/h);
gl.disable(gl.DEPTH_TEST);gl.disable(gl.BLEND);gl.useProgram(this.volume);gl.bindBuffer(gl.ARRAY_BUFFER,this.quad);let p=gl.getAttribLocation(this.volume,'position');gl.enableVertexAttribArray(p);gl.vertexAttribPointer(p,2,gl.FLOAT,false,0,0);
const u=n=>gl.getUniformLocation(this.volume,n);for(const [n,v] of [['eye',c.eye],['right',c.r],['up',c.u],['forward',c.f]])gl.uniform3fv(u(n),v);
gl.uniform1f(u('aspect'),w/h);gl.uniform1f(u('tanHalf'),c.t);gl.uniform1f(u('jmax'),this.data.jmax);gl.uniform1f(u('density'),this.mode===0?state.density:.07);gl.uniform2fv(u('extent'),this.data.extent);gl.uniform1i(u('cutaway'),this.mode===0&&state.cutaway?1:0);
gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,this.texture);gl.uniform1i(u('field'),0);gl.drawArrays(gl.TRIANGLES,0,6);gl.disableVertexAttribArray(p);
gl.useProgram(this.line);gl.bindBuffer(gl.ARRAY_BUFFER,this.buffer);for(const [name,size,offset] of [['position',3,0],['normal',3,12],['value',1,24]]){let a=gl.getAttribLocation(this.line,name);gl.enableVertexAttribArray(a);gl.vertexAttribPointer(a,size,gl.FLOAT,false,28,offset);}
const l=n=>gl.getUniformLocation(this.line,n);gl.uniformMatrix4fv(l('camera'),false,c.matrix);gl.uniform1i(l('mode'),this.mode);gl.uniform1i(l('cutaway'),state.cutaway?1:0);gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE);
gl.uniform1f(l('width'),this.mode===0?.007:.009);gl.uniform1f(l('opacity'),this.mode===0?.06:.06);gl.drawArrays(gl.TRIANGLES,0,this.count);
gl.enable(gl.DEPTH_TEST);gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);gl.uniform1f(l('width'),this.mode===0?.0012:.0028);gl.uniform1f(l('opacity'),this.mode===0?.5:.95);gl.drawArrays(gl.TRIANGLES,0,this.count);
for(const name of ['position','normal','value'])gl.disableVertexAttribArray(gl.getAttribLocation(this.line,name));gl.disable(gl.BLEND);}
}
