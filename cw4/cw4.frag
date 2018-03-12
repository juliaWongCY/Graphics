#version 400

in vec3 dir; 
out vec4 outcolour;

uniform mat4 mMatrix;
uniform mat4 mvMatrix;
uniform mat4 mvMatrixScene;
uniform mat4 pMatrix;
uniform mat3 normalMatrix; //mv matrix without translation

#define M_PI 3.1415926535897932384626433832795


uniform float shininess = 7;
uniform float ambientCoefficient = 0.3;
uniform float diffuseCoefficient = 1.3;
uniform float specularCoefficient = 1;
uniform float phi = 1000;

uniform float division_width = 0.4;
uniform float epsilon = 1e-4;
const vec3 light_source = vec3(6, 4, 3);
const float fresnel = 0.8;// Base on Schlick's approximation 

const int raytraceDepth = 42;
const int numSpheres = 6;

//example data structures
struct Ray
{
	vec3 origin;
	vec3 dir;
};

struct Sphere
{
	vec3 centre;
	float radius;
	vec3 colour;
};
struct Plane
{
	vec3 point;
	vec3 normal;
	vec3 colour;
};

struct Intersection
{
    float t; //closest hit
    vec3 point;	// hit point
    vec3 normal;	// normal
    int hit;	//did it hit?
    vec3 colour; // colour accumulation, can be also implemented in struct Ray
};

bool solveQuadratic(float a, float b, float c, inout float mui_0, inout float mui_1)
{
	// Finding the discriminant
	float dis =(b*b) - (4*a*c);
	if (dis < 0) return false;

	mui_0 = (-b - sqrt(dis)) / (2*a);
	mui_1 = (-b + sqrt(dis)) / (2*a);
		
	return true;
}

vec3 reflect(vec3 dir, vec3 normal) {
	//v' = v - (2v . n)n
	return dir - (dot(2*dir,normal) * normal);
}

vec3 refract(vec3 dir, vec3 normal, float ior) {
	//ior = index of refraction 
	// Need to make sure you know whether you entering/ leaving the transparent medium
	vec3 n = normal;
	float n_dot_v = dot(dir, n);
	float ior_1 = 1; //when the ray is entering the medium 
	float ior_2 = ior;
	
	if (n_dot_v < 0)  {
	  // The ray is entering the surface
	  n_dot_v = - n_dot_v;
	} else {
	  //The ray is already in the surface and leaving the surface
	  float temp = ior_1;
	  ior_1 = ior_2;
	  ior_2 = temp;
	  n = -n;
	}
	float coef = ior_1 / ior_2;
	float n_dot_v_sq = pow(n_dot_v, 2);
	float coef_sq = pow(coef,2);
	
	if ( n_dot_v_sq < (1 - coef_sq)) {
	  return vec3(0,0,0);
	} else {
	  float sq = sqrt(n_dot_v_sq+coef_sq-1) - n_dot_v;
	  return coef * (sq * n + dir);
	}
}

void computeBlinnPhongShading(Ray ray, inout Intersection intersect){
	float d = length(light_source - intersect.point);
	float c = phi / (4 * M_PI * d * d);
	vec3 n = normalize(intersect.normal);
	
	//Diffuse
	vec3 l = normalize(light_source - intersect.point);
	float intensity = max(dot(n, l), 0.0);
	vec3 diffuse_term = intersect.colour * diffuseCoefficient * intensity;
	
	//Specular
	//vec3 r = 2 * intensity * n - l;
	vec3 v = normalize(-ray.dir);
	vec3 h = normalize(1 + v);
	float r_c = pow(max(dot(n, h), 0.0), shininess);
	vec3 spec_term = intersect.colour * specularCoefficient * r_c;

	float shadow_factor = intersect.hit;

	intersect.colour = intersect.colour*ambientCoefficient + shadow_factor * (diffuse_term + spec_term) * c; 

	
}

void shpere_intersect(Sphere sph, Ray ray, inout Intersection intersect)
{
	float mu_0, mu_1;
	vec3 delta_p = ray.origin- sph.centre;
	// Preparing the quadratic formula, trying to solve mui
	// mui^2 + 2(mui)(dot(d, delta_p)) + (|delta_p ^ 2| - r^2) = 0
	// Define the coefficients a,b and c
	vec3 d = normalize(ray.dir); 
	float a = dot(d, d);
	float b = 2 * dot(d, delta_p);
	float c = dot(delta_p, delta_p) - pow(sph.radius, 2);
	if (solveQuadratic(a, b, c, mu_0, mu_1)) {

	  if (mu_0 < 0 && mu_1 < 0) {
		return;
	  } else if (mu_0 < 0) {
		intersect.t = mu_1;
	  } else if (mu_1 < 0) {
		intersect.t = mu_0;
	  } else{
	    intersect.t = min(mu_0, mu_1);
	  }
	    intersect.hit = 1;
	    intersect.point = ray.origin + intersect.t * ray.dir;
	    intersect.normal = normalize(intersect.point - sph.centre);
	    intersect.colour = sph.colour;
	  
	}
}

void plane_intersect(Plane pl, Ray ray, inout Intersection intersect)
{
	vec3 d = normalize(ray.dir);
	vec3 n = pl.normal;
	float denom = dot(d, n);
	float mu = - dot((ray.origin - pl.point), n) / denom;
	if (mu >= 0) {
  	  intersect.hit = 1;
	  intersect.t = mu;
	  intersect.normal = pl.normal;
	  intersect.point = ray.origin + d * intersect.t; // p_0 + d(mu)

	  bool x = int((intersect.point.x) / division_width) % 2 == 0;
	  //bool y = int((intersect.point.y) / division_width) % 2 == 0;
	  bool z = int((intersect.point.z) / division_width) % 2 == 0;
	
	  if (x ^^ z) {
		intersect.colour = pl.colour - pl.colour; 
	  } else {
	    intersect.colour = pl.colour;
	  }

	}	
}

Sphere sphere[numSpheres];
Plane plane;
void Intersect(Ray r, inout Intersection i)
{
	plane_intersect(plane, r, i);
	for (int s = 0; s < numSpheres; s++){
	  	shpere_intersect(sphere[s], r, i);
	}

/*
	Intersection temp;
	for (int s = 0; s < numSpheres; s++){
	  	shpere_intersect(sphere[s], r, temp);
	}
	plane_intersect(plane, r, i);
	float dist_s = distance(r.origin, temp.point);
	float dist_p = distance(r.origin, i.point);
	if ( dist_s < dist_p ) {
	  i = temp;
	}
	*/
}

int seed = 0;
float rnd()
{
	seed = int(mod(float(seed)*1364.0+626.0, 509.0));
	return float(seed)/509.0;
}

vec3 computeShadow(in Intersection intersect)
{
	Ray shadow;
	shadow.origin = intersect.point + epsilon * intersect.normal;
	shadow.dir = normalize(light_source - shadow.origin);
	// Compute shdow if there's some obstruction in between
	intersect.hit = 0;
	//Check if the shadow(secodary) ray is blocked by other object 
	Intersect(shadow, intersect);
	// When there's obstruction, then return 0 brightness
	if (intersect.hit == 1) {
	  return vec3(0,0,0);
	}
	return intersect.colour;
}

void main()
{
	
	//please leave the scene config unaltered for marking 
	sphere[0].centre   = vec3(-2.0, 1.5, -3.5);
	sphere[0].radius   = 1.5;
	sphere[0].colour = vec3(0.8,0.8,0.8);
	sphere[1].centre   = vec3(-0.5, 0.0, -2.0);
	sphere[1].radius   = 0.6;
	sphere[1].colour = vec3(0.3,0.8,0.3);
	sphere[2].centre   = vec3(1.0, 0.7, -2.2);
	sphere[2].radius   = 0.8;
	sphere[2].colour = vec3(0.3,0.8,0.8);
	sphere[3].centre   = vec3(0.7, -0.3, -1.2);
	sphere[3].radius   = 0.2;
	sphere[3].colour = vec3(0.8,0.8,0.3);
	sphere[4].centre   = vec3(-0.7, -0.3, -1.2);
	sphere[4].radius   = 0.2;
	sphere[4].colour = vec3(0.8,0.3,0.3);
	sphere[5].centre   = vec3(0.2, -0.2, -1.2);
	sphere[5].radius   = 0.3;
	sphere[5].colour = vec3(0.8,0.3,0.8);
	plane.point = vec3(0,-0.5, 0);
	plane.normal = vec3(0, 1.0, 0);
	plane.colour = vec3(1, 1, 1);
	seed = int(mod(dir.x * dir.y * 39786038.0, 65536.0));
	//scene definition end

	vec4 colour = vec4(0,0,0,1);

	Ray r;
	vec3 o = vec3(0,0,0);
	//As from the translation matrix (view matrix) --
	//The last column indicates the translation along the x, y and z axis respectively
	//(1 0 0 t_x)
	//(0 1 0 t_y)
	//(0 0 1 t_z)
	//(0 0 0  1	)
	float x = mvMatrixScene[3][0]; //left and right 
	float y = mvMatrixScene[3][1]; // up/down
	float z = mvMatrixScene[3][2] + 40.0 ; //zoom in/out
	r.origin = mat3(mvMatrixScene) * (o + vec3(x, y, z));
	r.dir = normalize(mat3(mvMatrixScene) * dir);

	for (int rd=0; rd < raytraceDepth; rd++) {
	  float contr = 1;
	  Intersection i;
	  i.hit = 0;
	  Intersect(r, i);
	  if (i.hit == 1) {
		//outcolour += i.colour;
	    contr = contr / (rd + 1);
		computeBlinnPhongShading(r, i);
		outcolour += contr * vec4(computeShadow(i),1);
	  } else {
	 	break;
	  }
	 
	  r.origin = i.point;
	  r.origin += epsilon * i.normal;
//	  r.dir = reflect(r.dir, i.normal);
	  // assume the balls are glasses
	  r.dir = fresnel * reflect(r.dir, i.normal) + (1 - fresnel) * refract(r.dir, i.normal, 1.3 ); 
	}
	outcolour = outcolour / 2 ;
}





























