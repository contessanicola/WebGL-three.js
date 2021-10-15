import * as THREE from './three.js-master/src/Three.js';
//import { OrbitControls } from 'https://cdn.skypack.dev/three@0.130.0/examples/jsm/controls/OrbitControls.js';
import { FlyControls } from './three.js-master/examples/jsm/controls/FlyControls.js';

export default class Draw {
	constructor(){
		this.renderer = new THREE.WebGLRenderer();
		this.renderer.setSize(window.innerWidth,window.innerHeight);
		document.getElementById('container').appendChild(this.renderer.domElement);
		
		this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight,1,100)

		this.camera.position.x = 0;
		this.camera.position.y = 0;
		this.camera.position.z = 3;

		this.scene = new THREE.Scene();
		this.scene.background = new THREE.Color(0xFFFF00);
		this.clock = new THREE.Clock();
		
		this.addQuad();
		this.time = 0;

		this.controls = new FlyControls(this.camera, this.renderer.domElement);
		this.controls.movementSpeed = 1;
		this.controls.domElement = this.renderer.domElement;
		this.controls.rollSpeed = 0.40;
		this.controls.autoForward = false;
		this.controls.dragToLook = true;
		this.test = 1;
		document.addEventListener("keydown", this.onDocumentKeyDown.bind(this), false);
		//this.resize();
		this.render();
		//this.setupResize();
	}

	addQuad(){
		this.geometry = new THREE.PlaneGeometry(2,2);
		this.shaderMaterial = new THREE.ShaderMaterial({
			vertexShader:  
			`
			void main() {
				//gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); 
				gl_Position = vec4(position, 1.0); 
			}
			`,
			fragmentShader: 
			`
			uniform vec2 iResolution;
			uniform float iTime;
			uniform int iScene;
			// = object.matrixWorld
			uniform mat4 modelMatrix;
			// = camera.matrixWorldInverse * object.matrixWorld
			uniform mat4 modelViewMatrix;
			// = camera.projectionMatrix
			uniform mat4 projectionMatrix;
			// = inverse transpose of modelViewMatrix
			uniform mat3 normalMatrix;

			#define BACKGROUND_COLOR vec3(0.5,0.8,0.9)//vec3(0.6,0.8,1.0)
			#define FOG_COLOR vec3(0.30, 0.36, 0.40)
			#define LIGHT_COLOR vec3(1.0,0.9,0.7)
			#define LIGHT_DIRECTION vec3(0.36, 0.48, 0.80)

			#define MAX_DIST 500.
			#define MIN_DIST 1e-4
			#define MAX_MARCHES 500
			#define SUN_SIZE 0.001
			#define SUN_SHARPNESS 1.5
			#define POWER 8.
			#define OSCILLATION 1
			#define AA 1

			vec4 mandelbulbTrap;
			float res;
			
			float sdSphere(vec3 p, float s){
				return length(p)-s;
			}
			
			float sdBox(vec3 p, vec3 b){
				vec3 q = abs(p) - b;
				return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
			  }

			float sdPlane(vec3 p){
				return p.y;
			} 	
			
			vec3 rotateY(vec3 p, float alpha){
				float px=p.x;
				float c=cos(alpha);
				float s=sin(alpha);
				
				 p.x=c*px-s*p.z;
				p.z=s*px+c*p.z;
				
				return p;
			}
			
			vec3 rotateX(vec3 p, float alpha){
				float py=p.y;
				float c=cos(alpha);
				float s=sin(alpha);
				
				 p.y=c*py-s*p.z;
				p.z=s*py+c*p.z;
				
				return p;
			}

			float sdMenger(vec3 p){ //https://www.shadertoy.com/view/wtfcWN by flo72
				float size=2.;
				p.z -=3.;  
				//p=rotateY(p,iTime*.5);
				vec3[] s = vec3[](vec3(1,1,1),vec3(1,1,0));
				
				for(int iter=0;iter<5;++iter){
					float alpha=(sin(iTime*0.1)+1.)  ; 
					p=rotateY(p,alpha);
					//float beta=(sin(iTime*0.1)+1.)  ; 
					//p=rotateX(p,beta);
				   
					p=abs(p);
					if(p.y > p.x) p.yx = p.xy;
					if(p.z > p.y) p.zy = p.yz;
					
					if(p.z > .5*size) p -= size*s[0];
					else p -= size*s[1];
					size /=3.;
					
				}
				return sdBox(p,vec3(1.5*size));
			}
			
			float sdMandelbulb(vec3 p) {
				vec3 z = p;
				
				float m = dot(z,z);
				vec4 trap = vec4(abs(z),m);
				float power;
				if(OSCILLATION == 1)
					power = (sin(iTime*0.1)+1.) * 8./2. + 1.;
				float dr = 1.0;
				float r = 0.0;
				for (int i = 0; i < 15 ; i++) {
					r = length(z);
					if (r>2.) break;	
					// convert to polar coordinates
					float theta = power*acos(z.z/r);
					float phi = power*atan(z.y,z.x);
					dr =  power*pow( r, power-1.)*dr + 1.0;	
					// scale and rotate the point	
					// convert back to cartesian coordinates
					z = pow(r,power)*vec3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
					//w = p + pow(r,8.0) * vec3( sin(b)*sin(a), cos(b), sin(b)*cos(a) );
					
					z+=p;
			
					trap = min(trap, vec4(abs(z),m) );
					m = dot(z,z);					
				}
				mandelbulbTrap = vec4(m,trap.yzw);
				return 0.5*log(r)*r/dr;
			}

			float sdSphereMod(vec3 p, float s){
				vec3 sphere = vec3 (1.0,1.0,1.0);
				return length(mod(sphere.xyz -p,s) - vec3(s/2.0)) - .5;
				return length(p)-s;
			}

			float distanceField(vec3 p){
				if(iScene == 1){
					float Menger = sdMenger(p-vec3(0.0,0.0,-15.));
					return Menger;					
				}
				else if(iScene == 2){
					float Box = sdBox(p,vec3(1.5));
					if(Box > 0.1) return Box;
					float Mandelbulb = sdMandelbulb(p);
					return Mandelbulb;				
				}
				else if(iScene == 3){
					float SphereMod = sdSphereMod(p,2.0);
        			return SphereMod;
					
				}
				else if(iScene == 4){
					float Sphere = sdSphere(p-vec3(0.0,0.0,0.0),2.0);
					float Plane = sdPlane(p-vec3(0.0,-2.0,0.0));
					return min(Sphere,Plane);
				}
				else if(iScene == 5){
				}

				/**/

				
			}

			vec3 calcNormal(vec3 p, float h){ // https://www.iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
				vec3 k = vec3(1,-1,0);
				return normalize( k.xyy*distanceField( p + k.xyy *h) + 
								  k.yyx*distanceField( p + k.yyx *h) + 
								  k.yxy*distanceField( p + k.yxy *h) + 
								  k.xxx*distanceField( p + k.xxx *h) );
			}
			/*
			float softShadow(in vec3 ro, in vec3 rd, float mint, float maxt, float k, float min_dist){ //https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
				float res_sha = 1.0;
				float ph = 1e10;
				for(float t=mint; t<maxt;){
					float h = distanceField(ro + rd*t);
					if( h < min_dist)
						return 0.0;
					float y = h*h/(2.0*ph);
					float d = sqrt(h*h-y*y);
					res_sha = min(res_sha, k*d/max(0.0,t-y));
					ph = h;
					t += h;  
				}
				res_sha = clamp( res_sha, 0.0, 1.0 );
				return res_sha*res_sha*(3.0-2.0*res_sha);
				//return res_sha;
			}*/

			float shadow( in vec3 ro, in vec3 rd, float mint, float maxt, float k, float min_dist)
			{
				float res_sha = 1.0;
				for( float t=mint; t<maxt; )
				{
					float h = distanceField(ro + rd*t);
					if( h<min_dist )
						return 0.0;
						res_sha = min( res_sha, k*h/t );
					t += h;
				}
				res_sha = clamp( res_sha, 0.0, 1.0 );
				return res_sha*res_sha*(3.0-2.0*res_sha);
			}
			
			float ambientOcclusion(vec3 p, vec3 n){
				float steps = 0.1;
				float ao = 0.0;
				float dist;
				for(float i=1.; i<= 3.;i+=1.){
					dist = steps * i;
					ao += max(0.0,(dist - distanceField(p+n*dist))/ dist);
				}
				return (1.0-ao * 0.22);
			}

			vec3 raymarching(vec3 ro, vec3 rd, out int iter){
				float min_d = 1.0;
				float t = 0.; //distance traveled alongside the ray vector
				float s = 0.;
				float d = 0.;
				for(int i = 0; i < MAX_MARCHES; i++){   
					iter = i;    
					float min_dist = max(res*t, MIN_DIST);
					//float min_d = min(min_d, 10.0 * d / t); Check how close we got without hitting so we can use it for something eg GLOW EFFECT
					d = distanceField(ro + rd * t);
			
					if (t > MAX_DIST){ break; }
					else if (d < min_dist){
						s += d / min_dist;         // can use this for ambient occlusion
						break;
					}     
					t += d;      
				}
				return vec3(d,t,s);
			} 

			vec4 render(vec3 ro, vec3 rd){
				vec4 col = vec4(0.0);
				int iter;
				vec3 raymarch = raymarching(ro,rd,iter);
				float d = raymarch.x;
				float t = raymarch.y;
				float s = raymarch.z;
			
				float min_dist = max(res*t, MIN_DIST);
				vec3 p = ro + rd * t;
				if(d < min_dist){
					vec3 n = calcNormal(p,min_dist);
					float ks = 1.0;
					if(iScene == 1)
						col.xyz = vec3(0.2125, 0.1875, 0.09);
					else
						col.xyz = vec3(0.01);
					col.xyz = mix( col.xyz, vec3(0.10,0.20,0.30), clamp(mandelbulbTrap.y,0.0,1.0) );
					col.xyz = mix( col.xyz, vec3(0.02,0.10,0.30), clamp(mandelbulbTrap.z*mandelbulbTrap.z,0.0,1.0) );
					col.xyz = mix( col.xyz, vec3(0.30,0.10,0.02), clamp(pow(mandelbulbTrap.w,6.0),0.0,1.0) );
					col.xyz *= 0.5;

					vec3 sun_light = LIGHT_COLOR * clamp(dot(n,LIGHT_DIRECTION ), 0., 1.);//LIGHT_COLOR * max(dot(n, LIGHT_DIRECTION), 0.0);
					vec3 sky_light = (BACKGROUND_COLOR*0.10)* clamp(0.5+0.5*dot(n,vec3(0.,1.,0.)), 0., 1.);
					vec3 bounce_light = (vec3(.06,.063,.07))* clamp(0.5+0.5*dot(n,vec3(0.,-1.,0.)), 0., 1.);


					vec3 sun_shadow;
					if(iScene == 2 || iScene == 3){
						sun_shadow = vec3(1.);
					}
					else{
						sun_shadow = LIGHT_COLOR * shadow(p, LIGHT_DIRECTION, 0.1, 30.0, 12.0, min_dist);
					}
					
					
					//vec3 ref = reflect(rd, n);
					vec3 sun_half = normalize(LIGHT_DIRECTION-rd);
					vec3 sun_specular = ks*pow(clamp(dot(n,sun_half),0.0,1.0),8.0)*sun_light*(0.04+0.96*pow(clamp(1.0+dot(sun_half,rd),0.0,1.0),5.0));

					float ao = ambientOcclusion(p,n); 

					
					
					vec3 sum = vec3(0.0);
					sum += 8. * sun_light * sun_shadow;  
					sum += 2. * sky_light * ao;     
					sum += 4.0 * bounce_light * ao;    
					col.xyz *= sum;

					col.xyz += sun_specular * LIGHT_COLOR  * 8. * sun_shadow;
				}
				else {
					vec3 sky = BACKGROUND_COLOR - max(rd.y,0.0)*0.5;
					float sun = dot(rd, LIGHT_DIRECTION) - 1.0 + SUN_SIZE;
					sun = min(exp(sun * SUN_SHARPNESS / SUN_SIZE), 1.0);
					col.xyz += sky;
					col.xyz += LIGHT_COLOR * sun;
				}
				//return vec4(vec3(float(iter)/128.0), 1.0); //check ray iteration
				return col;
			} 

			void main() {
				res = 1.0 / 1080.0;
				vec4 col = vec4(1.0);
				vec2 uv = (gl_FragCoord.xy+0.5*(-iResolution.xy))/iResolution.y;
				vec3 ro = vec3(0.0, 0.0, -8.0);
    			vec3 rd = vec3(uv, -1.0);

				vec4 far_4 = inverse(modelViewMatrix) * vec4(normalize(vec3(uv.x,uv.y,-1.47)), 1.0);
           	 	vec3 far_3 = far_4.xyz/far_4.w;
				rd = normalize(far_3.xyz - cameraPosition);


				col = render(cameraPosition,rd);
				col = pow(col, vec4(0.4545));

				gl_FragColor = vec4(col);

				//gl_FragColor = vec4(uv,0.0,1.0);
			}
			`,
			uniforms:
			{
				iResolution: {value: new THREE.Vector2()},
				iTime: {value: 1.0},
				iScene: {value: 1}
			}
			
		})
		
		this.shaderMesh = new THREE.Mesh(this.geometry,this.shaderMaterial);
		this.shaderMesh.frustumCulled = false;
		this.shaderMaterial.uniforms.iResolution.value.x = window.innerWidth;
		this.shaderMaterial.uniforms.iResolution.value.y = window.innerHeight;
		this.scene.add(this.shaderMesh);
	}
	
	onDocumentKeyDown(event) {
		var keyCode = event.which;
		if (keyCode == 49) {
			this.shaderMaterial.uniforms.iScene.value = 1;	
		} else if (keyCode == 50) {
			this.shaderMaterial.uniforms.iScene.value = 2;
		} else if (keyCode == 51) {
			this.shaderMaterial.uniforms.iScene.value = 3;
		} else if (keyCode == 52) {
			this.shaderMaterial.uniforms.iScene.value = 4;
		} else if (keyCode == 53) {
			this.shaderMaterial.uniforms.iScene.value = 5;
		}
	}
/*
	setupResize(){
		window.addEventListener("resize", this.resize.bind(this));
	}

	resize(){
		this.width = this.container.offsetWidth;
		this.height = this.container.offsetHeight;
		this.renderer.setSize(this.width, whis.height);
		this.camera.aspect = this.width / this.height;
	}*/

	render(){
		this.time++;
		const delta = this.clock.getDelta();
		//this.mesh.rotation.x += 0.01;
		//console.log(this.test);
		this.renderer.render(this.scene,this.camera);
		this.shaderMaterial.uniforms.iTime.value = this.clock.elapsedTime;
		window.requestAnimationFrame(this.render.bind(this));

		this.controls.movementSpeed = 20 * delta;
		this.controls.update( delta );
	}

	
}

new Draw();

