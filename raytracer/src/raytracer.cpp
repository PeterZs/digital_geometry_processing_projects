#include <vector>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <map>
#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <memory>

#include "raytracer.hpp"
#include "image.hpp"
#include "cstdio"
#include "kdtree.hpp"

// 2016 Version

static void write_vtk_mesh(const std::vector<Vector>& verts, const std::vector<int>& conn)
{
	FILE *fl;
	fl = fopen("mesh.vtk", "w"); assert(fl);
	const int n_tri = conn.size() / 3;

	// write the header
	fprintf(fl, "# vtk DataFile Version 2.0\n");
	fprintf(fl, "Cartel output mesh\n");
	fprintf(fl, "ASCII\n");
	fprintf(fl, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fl, "\n");

	// write the vertices
	fprintf(fl, "POINTS %d float\n", (int)verts.size());
	for (auto it = verts.begin() ; it != verts.end() ; ++it)
	{
		fprintf(fl, "%e %e %e \n", (*it)[0], (*it)[1], (*it)[2]);
	}
	fprintf(fl, "\n");

	// write the faces
	fprintf(fl, "CELLS %d %d \n", n_tri, n_tri*4);
	for (int it = 0 ; it < n_tri ; it++)
	{
		fprintf(fl, "3 %d %d %d \n", conn[it*3+0], conn[it*3+1], conn[it*3+2]);
	}
	fprintf(fl, "\n");

	// write the face types
	fprintf(fl, "CELL_TYPES %d \n", n_tri);
	for (int f = 0 ; f < n_tri ; f++)
	{
		fprintf(fl, "5\n");
	}

	// close the file
	fclose(fl);
}

static void write_vtk_bbox(const std::vector<const BoundingBox*>& bbs, const Matrix* trans)
{
	FILE *fl;
	fl = fopen("bboxtri.vtk", "w"); assert(fl);

	// write the header
	fprintf(fl, "# vtk DataFile Version 2.0\n");
	fprintf(fl, "Cartel output mesh\n");
	fprintf(fl, "ASCII\n");
	fprintf(fl, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fl, "\n");

	// write the vertices
	fprintf(fl, "POINTS %d float\n", (int)bbs.size()*8);
	for (auto it = bbs.begin() ; it != bbs.end() ; ++it)
	{
		const BoundingBox &bbox = **it;

		Vector minp, maxp;

		if (trans)
		{
			minp = bbox.minp;
			maxp = bbox.maxp;
			minp[3] = 1.;
			maxp[3] = 1.;
			minp = (*trans)*minp;
			maxp = (*trans)*maxp;
		}
		else
		{
			minp =  bbox.minp;
			maxp = bbox.maxp;
		}

		fprintf(fl, "%e %e %e \n", minp[0], minp[1], minp[2]);
		fprintf(fl, "%e %e %e \n", minp[0], minp[1], maxp[2]);
		fprintf(fl, "%e %e %e \n", minp[0], maxp[1], maxp[2]);
		fprintf(fl, "%e %e %e \n", minp[0], maxp[1], minp[2]);
		fprintf(fl, "%e %e %e \n", maxp[0], minp[1], minp[2]);
		fprintf(fl, "%e %e %e \n", maxp[0], minp[1], maxp[2]);
		fprintf(fl, "%e %e %e \n", maxp[0], maxp[1], maxp[2]);
		fprintf(fl, "%e %e %e \n", maxp[0], maxp[1], minp[2]);
	}
	fprintf(fl, "\n");

	// write the faces
	fprintf(fl, "CELLS %d %d \n", (int)bbs.size(), (int)bbs.size()*9);
	for (int ib = 0 ; ib < (int)bbs.size() ; ib++)
	{
		fprintf(fl, "8 ");
		for (int i = 0 ; i < 8 ; i++)
			fprintf(fl, "%d ", ib*8+i);

		fprintf(fl, "\n ");
	}
	fprintf(fl, "\n");

	// write the face types
	fprintf(fl, "CELL_TYPES %d \n", (int)bbs.size());
	for (int f = 0 ; f < (int)bbs.size() ; f++)
	{
		fprintf(fl, "12\n");
	}

	// close the file
	fclose(fl);
}


void Raytracer::render(const char *filename, const char *depth_filename, Scene const &scene)
{
	/*
	 * Allocate the two images that will ultimately be saved.
	 */
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);

    /*
     * Create the zBuffer.
     */
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++)
    {
        zBuffer[i] = DBL_MAX;
    }

    /*
     * Calculate the camera properties
     * I will opt for a right handed coordinate system.
     *    Yscreen------------------------
     *    ^                             |
     *    |               ^             |
     *    |               | v_cam       |
     *    |               |             |
     *    |               |             |
     *    |      <---------             |
     *    |       u_cam                 |
     *    |                             |
     *    |                             |
     *    -----------------------------> Xscreen
     */

    const Vector &cam_loc = scene.camera.position;

    // principle directions
    const Vector cam_w = (scene.camera.center - cam_loc).normalized();
    const Vector cam_v = (scene.camera.up - scene.camera.up.dot(cam_w)*cam_w).normalized();
    const Vector cam_u = (cam_v.cross(cam_w)).normalized();


    // frustum
    const double cam_near = scene.camera.zNear;
    const double cam_far = scene.camera.zFar;
    const double cam_height = 2. * tan(deg2rad(scene.camera.fov/2.)) * cam_near;
    const double cam_width = scene.camera.aspect * cam_height;

    // left upper corner of the projection screen
    const Vector dw = cam_w * cam_near;
    const Vector du = (cam_width  / 2.) * cam_u;
    const Vector dv = (cam_height / 2.) * cam_v;

    const Vector screen_corner1 = cam_loc + dw + du + dv;
    const Vector screen_corner2 = cam_loc + dw - du + dv;
    const Vector screen_corner3 = cam_loc + dw - du - dv;
    const Vector screen_corner4 = cam_loc + dw + du - dv;

    /*
     * Debug
     * write the camera and the screen in a vtk file.
     */
    if(true)
    {
    	std::cerr << "Near: " << cam_near << std::endl;
      	std::cerr << "Far: " << cam_far << std::endl;
      	std::cerr << "Height: " << cam_height << std::endl;
      	std::cerr << "Width: " << cam_width << std::endl;
     	std::cerr << "Ray centered at: " << cam_loc << std::endl;
     	std::cerr << "Corner1: " << screen_corner1  << std::endl;
     	std::cerr << "Corner2: " << screen_corner2  << std::endl;
     	std::cerr << "Corner3: " << screen_corner3  << std::endl;
     	std::cerr << "Corner4: " << screen_corner4  << std::endl;

     	/*
     	 * Write camera and meshes
     	 */
    	std::vector<Vector> verts =
    	{
    	//		cam_loc,
		//		screen_corner1,
		//		screen_corner2,
		//		screen_corner3,
		//		screen_corner4,
    	};
    	std::vector<int> conn =
    	{
    	//		0, 3, 2,
		//		0, 4, 3,
		//		4, 0, 1,
		//		1, 0, 2,
		//		3, 1, 2,
		//		3, 4, 1
    	};

    	int nbeg = verts.size();
    	std::vector<const Mesh*> meshes;
    	for (auto it = scene.objects.begin() ; it != scene.objects.end() ; ++it)
    	{
    		const Mesh *mesh = dynamic_cast<const Mesh*>(*it);
    		if(!mesh) continue;

    		meshes.push_back(mesh);
    		for (auto ip = mesh->positions.begin() ; ip!=mesh->positions.end() ; ip++)
    		{
    			Vector loc = (*ip);
    			loc[3] = 1.;
    			verts.push_back(mesh->transform*loc);
    		}

    		for (uint tri_i = 0; tri_i < mesh->triangles.size(); tri_i++)
    		{
    			Triangle const &tri = mesh->triangles[tri_i];
    			conn.push_back(tri[0].pi+nbeg);
    			conn.push_back(tri[1].pi+nbeg);
    			conn.push_back(tri[2].pi+nbeg);
    		}

     		nbeg += mesh->positions.size();
    	}
    	write_vtk_mesh(verts, conn);

    	/*
    	 * Write the bounding box for meshes.
    	 */
    	if(meshes.size())
    	{
    		const Mesh* mesh = meshes.front();
    		BoundingBox meshbbox(mesh->bboxMin, mesh->bboxMax);
    		std::vector<const BoundingBox*> treebboxs;
    		mesh->treeroot->get_members_bbox(treebboxs);
    		write_vtk_bbox(treebboxs, &mesh->transform);
    	}

    	// don't raytrace for now
    //	delete[] zBuffer;
    //	return;
    }

    /*
     * Iterate over all the pixels in the image.
     */
    for(int y = 0; y < scene.resolution[1]; y++)
    {
        for(int x = 0; x < scene.resolution[0]; x++)
        {

            /*
             *  Generate the appropriate ray for this pixel
             */
			Ray ray;
			if (scene.objects.empty())
			{
				//no objects in the scene, then we render the default scene:
				//in the default scene, we assume the view plane is at z = 640 with width and height both 640
				ray = Ray(scene.camera.position, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - scene.camera.position).normalized());
			}
			else
			{
				const double du = double(x+0.5) / scene.resolution[XX] * cam_width;
				const double dv = double(y+0.5) / scene.resolution[YY] * cam_height;
				const Vector pixel_cent = screen_corner4 - du * cam_u + dv * cam_v;
				ray = Ray(cam_loc, (pixel_cent-cam_loc).normalized());
			}

            // Initialize recursive ray depth.
            int ray_depth = 0;
           
            // Our recursive raytrace will compute the color and the z-depth
            Vector color;

            // This should be the maximum depth, corresponding to the far plane.
            // NOTE: This assumes the ray direction is unit-length and the
            // ray origin is at the camera position.
            double depth = scene.camera.zFar;

            // Calculate the pixel value by shooting the ray into the scene
            trace(ray, scene, &depth, &color, &ray_depth);

            // Depth test
            if((depth >= scene.camera.zNear) && (depth <= scene.camera.zFar) &&
               (depth < zBuffer[x + y*scene.resolution[0]]) )
            {
                zBuffer[x + y*scene.resolution[0]] = depth;

                // Set the image color (and depth)
                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) / 
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

		//output step information
		if (y % 100 == 0)
		{
			printf("Row %d pixels done.\n", y);
		}
    }

	//save image
    colorImage.writeBMP(filename);
    depthImage.writeBMP(depth_filename);

	printf("Ray tracing terminated and images are saved.\n");

    delete[] zBuffer;
}


bool Raytracer::trace(Ray const &ray,
		Scene const &scene,
		double *depth,
		Vector *rayOutColor,
		int *ray_depth)
{
    // Increment the ray depth.
	if(ray_depth) (*ray_depth)++;

	bool did_hit = false;

    // - iterate over all objects calling Object::intersect.
    // - don't accept intersections not closer than given depth.
    // - call Raytracer::shade with the closest intersection.
    // - return true iff the ray hits an object.
	if (scene.objects.empty())
	{
		// no objects in the scene, then we render the default scene:
		// For default, we assume there's a cube centered on (0, 0, 1280 + 160) with side length 320 facing right towards the camera
		// test intersection:
		double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
		double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
		if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
		{
			//if intersected:
			Material m; m.emission = Vector(16.0, 0, 0); m.reflect = 0; //just for default material, you should use the intersected object's material
			Intersection intersection;	//just for default, you should pass the intersection found by calling Object::intersect()
			*rayOutColor = shade(ray, *ray_depth, intersection, m, scene);
			intersection.depth = 1280;	//the depth should be set inside each Object::intersect()
		}
	}
	else
	{
		Intersection intersec(*depth);
		Material const *material = NULL;

		for (auto it = scene.objects.begin() ; it != scene.objects.end() ; ++it)
		{
			const Object *obj = *it;
			if(obj->intersect(ray, intersec))
			{
				did_hit = true;
				material = &obj->material;
			}
		}

		if(did_hit)
		{
			*depth = intersec.depth;
			if(ray_depth) *rayOutColor = shade(ray, *ray_depth, intersec, *material, scene);
		}
		else
		{
			*depth = intersec.depth;
			if(ray_depth) *rayOutColor = Vector(0.);
		}
	}

    // Decrement the ray depth.
	if(ray_depth) (*ray_depth)--;

    return did_hit;
}


Vector Raytracer::shade(Ray const &ray, int &ray_depth, Intersection const &intersection, Material const &material, Scene const &scene)
{
    // - iterate over all lights, calculating ambient/diffuse/specular contribution
    // - use shadow rays to determine shadows
    // - integrate the contributions of each light
    // - include emission of the surface material
    // - call Raytracer::trace for reflection/refraction colors
    // Don't reflect/refract if maximum ray recursion depth has been reached!

	/*
	 * First find the effect of the lights directly on the desired point.
	 */
	Vector diffuse(0), ambient(0), specular(0);
	const Vector &normal_vector = intersection.normal;
	const Vector view_vector = (ray.origin - intersection.position).normalized();

	for (auto lit = scene.lights.begin(); lit != scene.lights.end(); lit++)
	{
		const PointLight &light = *lit;

		/*
		 * Find the vectors to the light
		 */
		Vector light_vector = light.position - intersection.position;
		const double dist = light_vector.length();
		light_vector.normalize();
		Vector reflec_vector = -light_vector +
				2.0 * normal_vector.dot(light_vector) * normal_vector;
		reflec_vector.normalize();
		const double atten = light.attenuation[0]+
				light.attenuation[1]*dist+
				light.attenuation[2]*dist*dist;

		/*
		 * Find the contribution of the light to the surface if there were
		 * no shadows.
		 */
		Vector diffuse_contrib(0), ambient_contrib(0), specular_contrib(0);
		for(int i = 0 ; i < 3 ; i++)
		{
			ambient_contrib[i] += material.ambient[i] * light.ambient[i];
			diffuse_contrib[i] +=  material.diffuse[i] * light.diffuse[i]/atten *
					std::max(light_vector.dot(normal_vector), 0.) ;
			specular_contrib[i] += material.specular[i] * light.specular[i]/atten *
					pow(std::max( reflec_vector.dot(view_vector) , 0.) , material.shininess);
		}


		/*
		 * Using a shadow ray find out if the light is masked or not.
		 * If it is, adjust the contributions accordingly.
		 */
		const double epsilon = 1e-5; //self intersection!
		double ptolight_depth = dist;
		Ray ptolight_ray(intersection.position+epsilon*light_vector, light_vector);
		const bool is_light_masked = trace(ptolight_ray, scene, &ptolight_depth);

		if(is_light_masked)
		{
			diffuse_contrib = (1. - material.shadow) * diffuse_contrib;
			specular_contrib = (1. - material.shadow) * specular_contrib;
		}

		/*
		 * Add the contibution of this light
		 */
		ambient += ambient_contrib;
		diffuse += diffuse_contrib;
		specular += specular_contrib;

	} // End of for over lights

	/*
	 * Now find the reflected color on the same surface
	 */
	Vector reflec_color;
	if ( (ABS_FLOAT(material.reflect) > 1e-6) && (ray_depth < _max_ray_depth) )
	{
		// Find the reflection vector (reflection of view! not light).
		Vector reflec_vector = -view_vector +
				2.0 * normal_vector.dot(view_vector) * normal_vector;
		reflec_vector.normalize();

		//Emit a ray and find it's color
		const double epsilon = 1e-10; //self intersection!
		double reflec_depth = 1e6; //scene.camera.zFar // I don't really know! :D
		Ray reflec_ray(intersection.position+epsilon*reflec_vector, reflec_vector);

		trace(reflec_ray, scene, &reflec_depth,&reflec_color, &ray_depth);
	}

	return material.emission + ambient + diffuse + specular + material.reflect * reflec_color;
}
