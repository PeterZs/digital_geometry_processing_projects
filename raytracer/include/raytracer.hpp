#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "scene.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
using namespace std;


class Raytracer 
{
private:
	static int constexpr _max_ray_depth = 10;

public:
    // Render the given scene using raytracing.
    // Store the resulting image and depth map in the given files.
    static void render(const char *filename, const char *depth_filename, Scene const &scene);

private:
    // Shoot a ray into the scene and calculate the color/depth values.
    // Params:
    //   ray -- the ray being shot into the scene
    //   ray_depth -- the recursion depth of the current ray being traced
    //   scene -- the scene into which the ray is being shot
    //   rayOutColor -- set to the calculated color for the ray on exit
    //   depth -- set to the closest intersection z-depth on exit
    //      NOTE: The depth value passed in is treated as an upper bound.
    //            Objects further away than this depth should be ignored.
    // Returns true iff the ray hits an object in the scene.
    static bool trace(Ray const &ray, Scene const &scene, double *depth, Vector *rayOutColor = NULL, int *ray_depth = NULL);

    // Compute the shading at a given point and normal.
    // Uses the given material as well as the light sources 
    // and the other objects in the scene for shadows and reflections.
    // Params:
    //   ray -- the ray that intersected a surface
    //   ray_depth -- the recursion depth of the current ray being traced
    //   intersection -- information about the ray-surface intersection
    //   material -- the material of the object at the intersection point
    //   scene -- the scene being rendered; contains lights, other objects
    // Returns the calculated color.
	static Vector shade(Ray const &ray, int &ray_depth, Intersection const &intersection, Material const &material, Scene const &scene);
};


#endif
