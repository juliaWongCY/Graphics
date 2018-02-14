#version 410 
#define M_PI 3.1415926535897932384626433832795

in vec3 vertex_worldSpace;
in vec3 normal_worldSpace;
in vec2 textureCoordinate_input;


uniform mat4 mvMatrix;
uniform mat4 pMatrix;
uniform mat3 normalMatrix; //mv matrix without translation

uniform vec4 lightPosition_camSpace; //light Position in camera space

uniform vec4 ambient;
uniform vec4 diffuse;
uniform vec4 specular;
uniform float shininess;
uniform float ambientCoefficent;
uniform float diffuseCoefficent;
uniform float specularCoefficent;

uniform float phi = 5000;

out data
{
  vec4 position_camSpace;
  vec3 normal_camSpace;
  vec2 textureCoordinate;
  vec4 color;
}vertexInOut;

//Vertex shader compute the vectors per vertex
void main(void)
{
  //Put the vertex in the correct coordinate system by applying the model view matrix
  vec4 vertex_camSpace = mvMatrix*vec4(vertex_worldSpace,1.0f); 
  vertexInOut.position_camSpace = vertex_camSpace;
  
  //Apply the model-view transformation to the normal (only rotation, no translation)
  //Normals put in the camera space
  vertexInOut.normal_camSpace = normalize(normalMatrix*normal_worldSpace);
  
  //we need to make sure that the normals and texture coordinates
  //aren't optimised away, 
  //so we have to use them somehow.
  //Uniforms and array objects that are nor used for 
  //the final output(!) are  removed during 
  //glsl compilation regardless if you assign them. 
  vec4 workaround = 
		vec4((vertexInOut.normal_camSpace.x + textureCoordinate_input.x)*0.0001, 0, 0, 1);
  
  //forwarding pure red as RGBA color
  //Try to use the normals as RGB color or the texture coordiantes!

  //Gouraud shading
  //Finding the light vector by subtracting the object position from the light position
  float d = length(lightPosition_camSpace - vertexInOut.position_camSpace);
  vec3 lightVector = vec3(normalize(lightPosition_camSpace - vertexInOut.position_camSpace));
  float intensity = max(dot(vertexInOut.normal_camSpace, lightVector), 0.0);
  float coefficient = phi / (4 * M_PI * d * d);
  float radiance = (diffuseCoefficent * intensity) * coefficient;
  vertexInOut.color = vec4(1.0, 0.0, 0.0, 1.0) * radiance + ambient;
  
  //a negligible contribution from normals and texcoords is added 
  //to ensure these array objects are not optimsed away 
  vertexInOut.color += workaround;
  
  //Texture coordinate
  vertexInOut.textureCoordinate = textureCoordinate_input;
  
  gl_Position = pMatrix * vertex_camSpace;
}


