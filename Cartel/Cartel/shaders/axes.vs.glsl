#version 330

// Copyright (c) Russell Gillette
// December 2013

layout (location = 0) in vec3 Position;
layout (location = 1) in vec3 Normal;

out Data {
    vec3 color;
	vec2 texCoord;
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

uniform mat3 NormalMatrix;      // we keep a MV matrix without the translation component to apply to vectors
uniform mat4 MVP;               // ModelViewProjection Matrix

void main()
{
    // determine vertex color
    vec3 tnorm    = normalize( NormalMatrix * Normal );
    vec3 s        = normalize( vec3(Light0.LPosition) ); // incident vector
    outData.color = dot(s, tnorm) * Material0.Ka * 0.4 + Material0.Ka * 0.6;

    gl_Position = MVP * vec4(Position, 1.0);
}