#version 330

layout(triangles) in;
layout(triangle_strip, max_vertices=3) out;

in Data {
    vec3 color;
    float selected;
} inData[];

out wireData {
    vec3 f_color;
    vec3 dist_v;
    vec3 selected;
    noperspective vec3 dist;
} outData;


in gl_PerVertex {
    vec4 gl_Position;
    float gl_PointSize;
    float gl_ClipDistance[];
} gl_in[];

uniform           vec2 scale;

void main(void)
{
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;
    vec3 p2 = gl_in[2].gl_Position.xyz;

    // heavily modified from 'http://strattonbrazil.blogspot.ca/2011/09/single-pass-wireframe-rendering_10.html'
    vec2 q0 = scale * p0.xy/gl_in[0].gl_Position.w;
    vec2 q1 = scale * p1.xy/gl_in[1].gl_Position.w;
    vec2 q2 = scale * p2.xy/gl_in[2].gl_Position.w;
    vec2 v0 = q2 - q1;
    vec2 v1 = q2 - q0;
    vec2 v2 = q1 - q0;
    float area = abs(v1.x * v2.y - v1.y * v2.x) / 2;

    vec3 e = vec3(length(v0), length(v1), length(v2));

    float l_01 = length(p1-p0);
    float l_12 = length(p2-p1);
    float l_02 = length(p0-p2);

    outData.selected = vec3(inData[0].selected, inData[1].selected, inData[2].selected);
    outData.dist_v  = vec3(0, l_01, l_02);
    outData.dist    = vec3(area/e[0],0,0);
    outData.f_color = inData[0].color;
    gl_Position     = gl_in[0].gl_Position;
    EmitVertex();

    outData.selected = vec3(inData[0].selected, inData[1].selected, inData[2].selected);
    outData.dist_v  = vec3(l_01, 0, l_12);
    outData.dist    = vec3(0,area/e[1],0);
    outData.f_color = inData[1].color;
    gl_Position     = gl_in[1].gl_Position;
    EmitVertex();

    outData.selected = vec3(inData[0].selected, inData[1].selected, inData[2].selected);
    outData.dist_v  = vec3(l_02, l_12, 0);
    outData.dist    = vec3(0,0,area/e[2]);
    outData.f_color = inData[2].color;
    gl_Position     = gl_in[2].gl_Position;
    EmitVertex();

    EndPrimitive();
}
