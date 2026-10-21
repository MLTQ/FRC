function createVolumeRenderer(canvas) {
 const gl=canvas.getContext('webgl2',{alpha:true,antialias:false,premultipliedAlpha:false});
 if(!gl)throw new Error('WebGL2 is unavailable; showing current samples.');
 const vertex=`#version 300 es
 in vec2 position; void main(){gl_Position=vec4(position,0.,1.);}`;
 const fragment=`#version 300 es
 precision highp float;
 uniform sampler2D currents;
 uniform highp sampler3D currents3;
 uniform int mode;
 uniform vec3 lo3, hi3, dims3;
 uniform vec2 resolution, dimensions;
 uniform vec4 bounds;
 uniform vec3 rightAxis, upAxis, forwardAxis;
 uniform float pixelsPerUnit, maximum, opacity, threshold;
 uniform int component;
 out vec4 result;
 vec3 spectrum(float t){
   vec3 a=vec3(.06,.12,.55), b=vec3(.0,.70,.85), c=vec3(.25,.85,.38);
   vec3 d=vec3(1.,.88,.15), e=vec3(1.,.30,.06), f=vec3(.78,.025,.015);
   float x=clamp(t,0.,1.)*5.;
   if(x<1.)return mix(a,b,x);if(x<2.)return mix(b,c,x-1.);
   if(x<3.)return mix(c,d,x-2.);if(x<4.)return mix(d,e,x-3.);return mix(e,f,x-4.);
 }
 vec4 sampleCurrent3(vec3 p){
   vec3 uvw=(p-lo3)/(hi3-lo3);
   if(any(lessThan(uvw,vec3(0.)))||any(greaterThan(uvw,vec3(1.))))return vec4(0.);
   vec3 q=vec3(uvw.x,uvw.z,uvw.y)*(dims3-1.);ivec3 ijk=ivec3(floor(q));vec3 f=fract(q);
   ijk=clamp(ijk,ivec3(0),ivec3(dims3)-2);
   vec4 c000=texelFetch(currents3,ijk,0),c100=texelFetch(currents3,ijk+ivec3(1,0,0),0),c010=texelFetch(currents3,ijk+ivec3(0,1,0),0),c110=texelFetch(currents3,ijk+ivec3(1,1,0),0);
   vec4 c001=texelFetch(currents3,ijk+ivec3(0,0,1),0),c101=texelFetch(currents3,ijk+ivec3(1,0,1),0),c011=texelFetch(currents3,ijk+ivec3(0,1,1),0),c111=texelFetch(currents3,ijk+ivec3(1,1,1),0);
   if(min(min(min(c000.w,c100.w),min(c010.w,c110.w)),min(min(c001.w,c101.w),min(c011.w,c111.w)))<.5)return vec4(0.);
   return mix(mix(mix(c000,c100,f.x),mix(c010,c110,f.x),f.y),mix(mix(c001,c101,f.x),mix(c011,c111,f.x),f.y),f.z);
 }
 vec4 sampleCurrent(vec3 p){
   if(mode==1)return sampleCurrent3(p);
   vec2 uv=vec2((p.x-bounds.z)/(bounds.w-bounds.z),length(p.yz)/bounds.y);
   if(any(lessThan(uv,vec2(0.)))||any(greaterThan(uv,vec2(1.))))return vec4(0.);
   vec2 q=uv*(dimensions-1.);ivec2 ij=ivec2(floor(q));vec2 f=fract(q);
   ij=clamp(ij,ivec2(0),ivec2(dimensions)-2);
   vec4 a=texelFetch(currents,ij,0),b=texelFetch(currents,ij+ivec2(1,0),0);
   vec4 c=texelFetch(currents,ij+ivec2(0,1),0),d=texelFetch(currents,ij+ivec2(1,1),0);
   if(min(min(a.w,b.w),min(c.w,d.w))<.5)return vec4(0.);
   return mix(mix(a,b,f.x),mix(c,d,f.x),f.y);
 }
 void main(){
   vec2 screen=(gl_FragCoord.xy-.5*resolution)/pixelsPerUnit;
   vec3 origin=screen.x*rightAxis+screen.y*upAxis;
   vec3 lo=mode==1?lo3:vec3(bounds.z,-bounds.y,-bounds.y),hi=mode==1?hi3:vec3(bounds.w,bounds.y,bounds.y);
   float scale=mode==1?.5*(hi3.y-lo3.y):bounds.y;
   vec3 direction=forwardAxis;
   vec3 safe=sign(direction+vec3(1e-12))*max(abs(direction),vec3(1e-7));
   vec3 a=(lo-origin)/safe,b=(hi-origin)/safe;
   vec3 near3=min(a,b),far3=max(a,b);
   float near=max(max(near3.x,near3.y),near3.z),far=min(min(far3.x,far3.y),far3.z);
   if(far<=near||maximum<=0.){result=vec4(0.);return;}
   float ds=(far-near)/384.;vec4 accumulated=vec4(0.);
   for(int i=0;i<384;i++){
     vec3 p=origin+(near+(float(i)+.5)*ds)*direction;
     vec4 v=sampleCurrent(p);if(v.w<.5)continue;
     float value=component==0?length(v.xyz):v[component-1];
     float t=clamp(abs(value)/maximum,0.,1.);if(t<=threshold)continue;
     vec3 color=component==0?spectrum(t):(value<0.?mix(vec3(.94),vec3(.15,.38,.95),t):mix(vec3(.94),vec3(.98,.19,.09),t));
     float alpha=1.-exp(-opacity*pow(t,1.6)*ds/scale);
     accumulated.rgb+=(1.-accumulated.a)*color*alpha;
     accumulated.a+=(1.-accumulated.a)*alpha;
     if(accumulated.a>.995)break;
   }
   result=vec4(accumulated.rgb,1.);
 }`;
 function compile(type,source){const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw new Error(gl.getShaderInfoLog(s));return s;}
 const program=gl.createProgram();gl.attachShader(program,compile(gl.VERTEX_SHADER,vertex));gl.attachShader(program,compile(gl.FRAGMENT_SHADER,fragment));gl.linkProgram(program);
 if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(program));
 gl.useProgram(program);
 const buffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array([-1,-1,3,-1,-1,3]),gl.STATIC_DRAW);
 const attribute=gl.getAttribLocation(program,'position');gl.enableVertexAttribArray(attribute);gl.vertexAttribPointer(attribute,2,gl.FLOAT,false,0,0);
 const texture=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,texture);
 gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);
 gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);
 const uniforms={};for(const name of ['currents','resolution','dimensions','bounds','rightAxis','upAxis','forwardAxis','pixelsPerUnit','maximum','opacity','threshold','component','currents3','mode','lo3','hi3','dims3'])uniforms[name]=gl.getUniformLocation(program,name);
 // 3D runs: an RGBA16F 3D texture (Jr, Jtheta, Jz, support), x fastest in memory order (z, y, x) reversed below
 const texture3=gl.createTexture();gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_3D,texture3);
 for(const k of [gl.TEXTURE_MIN_FILTER,gl.TEXTURE_MAG_FILTER])gl.texParameteri(gl.TEXTURE_3D,k,gl.NEAREST);
 for(const k of [gl.TEXTURE_WRAP_S,gl.TEXTURE_WRAP_T,gl.TEXTURE_WRAP_R])gl.texParameteri(gl.TEXTURE_3D,k,gl.CLAMP_TO_EDGE);
 gl.activeTexture(gl.TEXTURE0);let previous3=null;
 let previous=null;
 return {draw(frame,options){
   const {yaw,pitch,zoom,maximum,component,opacity,threshold}=options;
   gl.viewport(0,0,canvas.width,canvas.height);gl.useProgram(program);
   const three=!!frame.volume3d_data;gl.uniform1i(uniforms.mode,three?1:0);gl.uniform1i(uniforms.currents3,1);
   if(three){const v=frame.volume3d,[nx,ny,nz]=v.shape;
     if(previous3!==frame.volume3d_data){gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_3D,texture3);gl.texImage3D(gl.TEXTURE_3D,0,gl.RGBA16F,nz,ny,nx,0,gl.RGBA,gl.HALF_FLOAT,frame.volume3d_data);gl.activeTexture(gl.TEXTURE0);previous3=frame.volume3d_data;}
     gl.uniform3f(uniforms.lo3,v.bounds[0],v.bounds[1],v.bounds[2]);gl.uniform3f(uniforms.hi3,v.bounds[3],v.bounds[4],v.bounds[5]);gl.uniform3f(uniforms.dims3,nz,ny,nx);}
   else if(previous!==frame.volume_data){gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA32F,...frame.volume.shape,0,gl.RGBA,gl.FLOAT,frame.volume_data);previous=frame.volume_data;}
   gl.uniform1i(uniforms.currents,0);gl.uniform2f(uniforms.resolution,canvas.width,canvas.height);
   if(frame.volume){gl.uniform2fv(uniforms.dimensions,frame.volume.shape);gl.uniform4fv(uniforms.bounds,frame.volume.bounds);}
   const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);
   gl.uniform3f(uniforms.rightAxis,cy,sy,0);gl.uniform3f(uniforms.upAxis,sp*sy,-sp*cy,cp);gl.uniform3f(uniforms.forwardAxis,-cp*sy,cp*cy,sp);
   gl.uniform1f(uniforms.pixelsPerUnit,zoom*Math.min(canvas.width,canvas.height)*.37);
   gl.uniform1f(uniforms.maximum,maximum/(three?(frame.volume3d.value_scale||1):1));gl.uniform1i(uniforms.component,component);
   gl.uniform1f(uniforms.opacity,opacity);gl.uniform1f(uniforms.threshold,threshold);gl.drawArrays(gl.TRIANGLES,0,3);
 }};
}
