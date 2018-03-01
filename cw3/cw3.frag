#version 410

uniform sampler2D textureRendered;

in vec2 varyingTextureCoordinate;

out vec4 fragColor;

float s[12] = float[](-0.10568,-0.07568,-0.042158,
-0.02458,-0.01987456,-0.0112458,
0.0112458,0.01987456,0.02458,
0.042158,0.07568,0.10568);

uniform float d_max = 0.3;
uniform int n = 12;
uniform vec2 c = vec2(0.5, 0.5);

void main(void)
{
  vec4 temp_acc;

  for (int i = 0; i <n ; i++) {
	float d_i = s[i] * d_max;
	vec2 p = varyingTextureCoordinate.st;
	vec2 p_vec = normalize(c - p);
    temp_acc = temp_acc + texture(textureRendered, (p + p_vec*d_i));
  }

  vec4 rgb_blur = temp_acc / n;	
  //Render the texture on a quad
  //fragColor = texture(textureRendered, varyingTextureCoordinate.st);
  fragColor = rgb_blur;

}




