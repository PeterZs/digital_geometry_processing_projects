#version 330

layout (location = 0) out vec4 FragColor;

in wireData {
    vec3 f_color;
    vec3 dist_v;
    vec3 selected;
    noperspective vec3 dist; // distances from edges
} inData;

// 0x1 -> faces
// 0x2 -> edges
// 0x4 -> verts
uniform int view_mode = 3; // default view mesh and edges

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
uniform MaterialInfo Material;

void main()
{

    // calculate edges
    float dst = min(min(inData.dist[0], inData.dist[1]), inData.dist[2]);
    float edgeIntensity = smoothstep(0, 1, dst);

    // calculate verts
    vec3 dist_v = inData.dist_v;
    int small = (dist_v[0] < dist_v[1]) ? 0 : 1;
    small = (dist_v[small] < dist_v[2]) ? small : 2;

    float threshold = 0.005; // change to make vertices appear bigger or smaller
    float vertIntensity = smoothstep(0, threshold, dist_v[small]);

    vec4 out_color = vec4(0, 0, 0, 0);

    if ((view_mode & 0x1) != 0) { // faces
        out_color = vec4(inData.f_color, 1.);
    }
    if ((view_mode & 0x2) != 0) { // edges
        out_color = mix(vec4(1, 1, 1, 1.), out_color, edgeIntensity);
    }
    if ((view_mode & 0x4) != 0) { // verts
        if (inData.selected[small] > 0.5)
            out_color = mix(vec4(1., 0., 0., 1.), out_color, vertIntensity);
        else
            out_color = mix(vec4(0., 0., 1., 1.), out_color, vertIntensity);
    }

    FragColor = out_color;
}