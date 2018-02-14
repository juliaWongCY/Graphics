#version 410

uniform vec4 ambient;
uniform vec4 diffuse;
uniform vec4 specular;
uniform float shininess;
uniform float ambientCoefficent;
uniform float diffuseCoefficent;
uniform float specularCoefficent;

uniform vec4 lightPosition_camSpace; //light Position in camera space

in fragmentData
{
  vec4 position_camSpace;
  vec3 normal_camSpace;
  vec2 textureCoordinate;
  vec4 color;
} frag;

out vec4 fragColor; 

//Fragment shader computes the final color
void main(void)
{
   //Toon shading
  vec3 l = vec3(normalize(lightPosition_camSpace - frag.position_camSpace));
  vec3 n = normalize(frag.normal_camSpace);
  float i_f = max(dot(n, l), 0.0);

  if (i_f > 0.98) {
    fragColor = vec4(0.8, 0.8, 0.8, 1.0);
  } else if (i_f > 0.5 && i_f <= 0.98){
	fragColor = vec4(0.8, 0.4, 0.4, 1.0);
  } else if (i_f > 0.25 && i_f <= 0.5) {
	fragColor = vec4(0.6, 0.2, 0.2, 1.0);
  } else {
	fragColor = vec4(0.1, 0.1, 0.1, 1.0);
  }

}


