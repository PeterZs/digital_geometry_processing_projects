#version 330

layout (location = 0) in vec3 Position;
layout (location = 1) in vec3 Normal;
layout (location = 2) in int  Selected;
layout (location = 3) in vec3 Color;

out Data {
	vec3 color;
    float selected;
} outData;

struct LightInfo
{
    vec4 LPosition; // we are going to treat this as a direction to achieve directional lighting
    vec3 La;       // ambient light
    vec3 Ld;       // diffuse light
    vec3 Ls;       // specular light
};
uniform LightInfo Light0;

struct MaterialInfo
{
    vec3 Ka;
    vec3 Kd;
    vec3 Ks;
    float Shininess;
};
uniform MaterialInfo Material0;

uniform mat4 ModelViewMatrix;
uniform mat3 NormalMatrix;      // we keep a MV matrix without the translation component to apply to vectors
uniform mat4 ProjectionMatrix;
uniform mat4 MVP;               // ModelViewProjection Matrix

void main()
{
    if (Color[0] >= 0) {
        outData.color = Color;
    } else {
        // determine vertex color
        vec3 tnorm     = normalize( NormalMatrix * Normal );
        vec3 s         = normalize( vec3(Light0.LPosition) ); // incident vector
		vec3 diffuse = dot(s, tnorm) * Material0.Kd;
		vec3 ambient = Material0.Ka;
		// we do not use a specular term
        outData.color = diffuse * 0.2 + ambient * 0.8;
    }

    if (Selected == 1) {
        outData.selected = 1.0;
	} else {
        outData.selected = 0.0;
	}
    gl_Position = MVP * vec4(Position, 1.0);
}