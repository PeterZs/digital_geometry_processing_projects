#version 330

layout (location = 0) out vec4 FragColor;

in Data {
    vec3 color;
	vec2 texCoord;
} inData;

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
    vec3 Ka;            // Ambient reflectivity
    vec3 Kd;            // Diffuse reflectivity
    vec3 Ks;            // Specular reflectivity
    float Shininess;    // Specular shininess factor
};
uniform MaterialInfo Material0;

void main()
{
    FragColor = vec4(inData.color,1.0);
}