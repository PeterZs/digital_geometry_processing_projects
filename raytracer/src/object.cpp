#include "object.hpp"
#include "kdtree.hpp"

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <cassert>


GetPot getpot;

// 2016 Version

bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assert the correct values of the W coords of the origin and direction.
    // You can comment this out if you take great care to construct the rays
    // properly.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray localRay(i_transform * ray.origin, i_transform * ray.direction);
	//!!! USEFUL NOTES: to calculate depth in localIntersect(), if the intersection happens at
	//ray.origin + ray.direction * t, then t is just the depth
	//!!! USEFUL NOTES: Here direction might be scaled, so you must not renormalize it in
	//localIntersect(), or you will get a depth in local coordinate system,
	//which can't be compared with intersections with other objects
    if (localIntersect(localRay, hit))
	{
        // Assert correct values of W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transform intersection coordinates into global coordinates.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}

//INTERSECTIONS 

bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
	const double a = ray.direction.dot(ray.direction);
	const double bp = ray.direction.dot(ray.origin);
	const double c = ray.origin.dot(ray.origin) - radius*radius;
	const double delta = bp*bp - a*c;

	// delta = 0 will almost never happen.
	// Have to take care of numerical artifacts using tolerance.
	if (delta < inter_tol) return false;

	const double sqrtdelta = sqrt(delta);
	const double t1 = (-bp + sqrtdelta) / a;
	const double t2 = (-bp - sqrtdelta) / a;

	double t;
	if (t1 > 0 && t2 > 0 ) t = std::min(t1,t2);
	else if (t1 > 0) t = t1;
	else  t = t2;

	if( (t > hit.depth) || (t < 0) ) return false;

	hit.normal = hit.position = ray.origin + t * ray.direction;
	hit.depth = t;
	assert( fabs(hit.position.dot(hit.position) - radius*radius) < 1e-10 );
	return true;
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
	if(fabs(ray.direction[ZZ]) < inter_tol) return false;

	const double t = -ray.origin[ZZ] / ray.direction[ZZ];

	if ((t > hit.depth) || (t < 0)) return false;

	// Find the normal
	hit.normal = Vector(0,0,1);
	if (ray.direction.dot(hit.normal, false) > 0)
		hit.normal = Vector(0,0,-1);

	// Find position and depth
	hit.position = ray.origin + ray.direction * t;
	hit.depth = t;

	return true;
}


Mesh::~Mesh()
{
	if(treeroot)
	{
		treeroot->free_children();
		delete treeroot;
		treeroot= NULL;
	}
}

bool Mesh::intersectTriangle(Ray const &ray, Triangle const &tri, Intersection &hit) const
{
	// Extract vertex positions from the mesh data.
	Vector const &p0 = positions[tri[0].pi];
	Vector const &p1 = positions[tri[1].pi];
	Vector const &p2 = positions[tri[2].pi];

	// Two of triangle edge
	const Vector e1 = p1 - p0;
	const Vector e2 = p2 - p0;

	// Find the equation of the surface
	const Vector vecareatot = e1.cross(e2);
	const Vector nsurf = vecareatot.normalized();

	// Find the intersection
	const double denom = ray.direction.dot(nsurf);
	if(fabs(denom) < 1e-13) return false;
	const double t = nsurf.dot(p0-ray.origin) / denom;

	// Find candidate point
	const Vector p = ray.origin + t * ray.direction;

	// Find if the ray is inside the triangle
	const Vector e3 = p0 - p;
	const Vector e4 = p1 - p;
	const Vector e5 = p2 - p;
	const Vector vecarea0 = e4.cross(e5);
	const Vector vecarea1 = e5.cross(e3);
	const Vector vecarea2 = e3.cross(e4);
	if (vecarea0.dot(vecareatot) < 0) return false;
	else if (vecarea1.dot(vecareatot) < 0) return false;
	else if (vecarea2.dot(vecareatot) < 0) return false;

	// alternative method - this is numerically unstable
	//if( (A0 + A1 + A2) > Atot ) return false;

	// Check depth
	if ((t > hit.depth) || (t < 0)) return false;

	// Interpolate the normal, direction, and depth
	const double A0 = vecarea0.length();
	const double A1 = vecarea1.length();
	const double A2 = vecarea2.length();
	const double Atot = vecareatot.length();
	const double alpha = A0/Atot;
	const double beta  = A1/Atot;
	const double gamma = A2/Atot;

	hit.depth = t;
	hit.normal = normals[tri[0].ni]*alpha + normals[tri[1].ni]*beta + normals[tri[2].ni]*gamma ;
	hit.position = p;

	return true;
}

bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
	//////////////////
	// YOUR CODE HERE (creative license)
    return false;
}


// Intersections!
// Use an octree for lightning speed-up + creative part.
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Hit the root of the tree
	return localIntersect_bbox(ray, treeroot, hit);
}

bool Mesh::localIntersect_bbox(Ray const &ray, TreeNode * subspace, Intersection &hit) const
{
	/*
	 * Check the subspace (bounding box)
	 */
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++)
	{
		// Ray parrallel to any of the bounding box faces
		if (ray.direction[i] == 0.0)
		{
			if (ray.origin[i] < subspace->bbox.minp[i] || ray.origin[i] > subspace->bbox.maxp[i])
			{
				return false;
			}
		}

		else
		{
			double t1 = (subspace->bbox.minp[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (subspace->bbox.maxp[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Ensure t1 <= t2

			if (t1 > tNear) tNear = t1; // We want the furthest tNear
			if (t2 < tFar) tFar = t2; // We want the closest tFar

			if (tNear > tFar) return false; // Ray misses the bounding box.
			if (tFar < 0) return false; // Bounding box is behind the ray.
		}
	}

	/*
	 * If the box has children recurse them, else check the members.
	 */
	if(subspace->left)
	{
		assert(subspace->right);
		const bool did_hit_right = localIntersect_bbox(ray, subspace->right, hit);
		const bool did_hit_left = localIntersect_bbox(ray, subspace->left, hit);
		return (did_hit_left || did_hit_right);
	}

	/*
	 * If there are no children check the member triangles
	 */
	if(!subspace->members.size()) return false;
	bool isHit = false;
	for (auto it = subspace->members.begin(); it != subspace->members.end(); ++it)
	{
		if (intersectTriangle(ray, **it, hit)) isHit = true;
	}
	return isHit;
}


double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}


bool Mesh::readOBJ(std::string const &filename)
{
	// Try to open the file.
	std::ifstream file(filename.c_str());
	if (!file.good()) {
		std::cerr << "Unable to open OBJ file \"" << filename << "\"" << std::endl;
		return false;
	}

	// Keep fetching op codes and processing them. We will assume that there
	// is one operation per line.
	while (file.good()) {

		std::string opString;
		std::getline(file, opString);

		std::stringstream opStream(opString);
		std::string opCode;
		opStream >> opCode;

		// Skip blank lines and comments
		if (!opCode.size() || opCode[0] == '#') {
			continue;
		}

		// Ignore groups.
		if (opCode[0] == 'g') {
			std::cerr << "ignored OBJ opCode '" << opCode << "'" << std::endl;

			// Vertex data.
		}
		else if (opCode[0] == 'v') {

			// Read in up to 4 doubles.
			Vector vec;
			for (int i = 0; opStream.good() && i < 4; i++) {
				opStream >> vec[i];
			}

			// Store this data in the right location.
			switch (opCode.size() > 1 ? opCode[1] : 'v') {
			case 'v':
				positions.push_back(vec);
				break;
			case 't':
				texCoords.push_back(vec);
				break;
			case 'n':
				normals.push_back(vec);
				break;
			case 'c':
				colors.push_back(vec);
				break;
			default:
				std::cerr << "unknown vertex type '" << opCode << "'" << std::endl;
				break;
			}

			// A polygon (or face).
		}
		else if (opCode == "f") {
			std::vector<Vertex> polygon;
			// Limit to 4 as we only can handle triangles and quads.
			for (int i = 0; opStream.good() && i < 4; i++) {

				// Retrieve a full vertex specification.
				std::string vertexString;
				opStream >> vertexString;

				if (!vertexString.size()) {
					break;
				}

				// Parse the vertex into a set of indices for position,
				// texCoord, normal, and colour, respectively.
				std::stringstream vertexStream(vertexString);
				std::vector<int> indices;
				for (int j = 0; vertexStream.good() && j < 4; j++) {
					// Skip slashes.
					if (vertexStream.peek() == '/') {
						vertexStream.ignore(1);
					}
					int index;
					if (vertexStream >> index)
						indices.push_back(index);
				}

				// Turn this into a real Vertex, and append it to the polygon.
				if (indices.size()) {
					indices.resize(4, 0);
					polygon.push_back(Vertex(
						indices[0] - 1,
						indices[1] - 1,
						indices[2] - 1,
						indices[3] - 1
						));
				}

			}

			// Only accept triangles...
			if (polygon.size() == 3) {
				triangles.push_back(Triangle(polygon[0],
					polygon[1],
					polygon[2]));
				// ...and quads...
			}
			else if (polygon.size() == 4) {
				// ...but break them into triangles.
				triangles.push_back(Triangle(polygon[0],
					polygon[1],
					polygon[2]));
				triangles.push_back(Triangle(polygon[0],
					polygon[2],
					polygon[3]));
			}

			// Any other opcodes get ignored.
		}
		else {
			std::cerr << "unknown opCode '" << opCode << "'" << std::endl;
		}
	}

	updateBBox();

	return true;
}

void Mesh::updateBBox()
{
	bboxMin = Vector(DBL_MAX, DBL_MAX, DBL_MAX);
	bboxMax = Vector(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	for (std::vector<Vector>::iterator pItr = positions.begin();
		pItr != positions.end(); ++pItr)
	{
		Vector const& p = *pItr;
		if (p[0] < bboxMin[0]) bboxMin[0] = p[0];
		if (p[0] > bboxMax[0]) bboxMax[0] = p[0];
		if (p[1] < bboxMin[1]) bboxMin[1] = p[1];
		if (p[1] > bboxMax[1]) bboxMax[1] = p[1];
		if (p[2] < bboxMin[2]) bboxMin[2] = p[2];
		if (p[2] > bboxMax[2]) bboxMax[2] = p[2];
	}

	// Update the bbox of each triangle as well
	for(auto it = triangles.begin() ; it != triangles.end() ; ++it)
	{
		it->bbox.minp = Vector(DBL_MAX, DBL_MAX, DBL_MAX);
		it->bbox.maxp = Vector(-DBL_MAX, -DBL_MAX, -DBL_MAX);
		for (int v = 0 ; v < 3 ; v++)
		{
			Vector &p = positions[(*it)[v].pi];
			if (p[0] < it->bbox.minp[0]) it->bbox.minp[0] = p[0];
			if (p[0] > it->bbox.maxp[0]) it->bbox.maxp[0] = p[0];
			if (p[1] < it->bbox.minp[1]) it->bbox.minp[1] = p[1];
			if (p[1] > it->bbox.maxp[1]) it->bbox.maxp[1] = p[1];
			if (p[2] < it->bbox.minp[2]) it->bbox.minp[2] = p[2];
			if (p[2] > it->bbox.maxp[2]) it->bbox.maxp[2] = p[2];
		}
	}
}

void Mesh::init()
{
	/*
	 * Set up kdtree (hierarchy of bounding boxes)
	 */
	int depth_max = getpot.follow(10, "-d");
	std::vector<const Triangle*> tripts;
	for (auto it = triangles.begin(); it != triangles.end() ; ++it)
	{
		tripts.push_back(&(*it));
	}
	treeroot = TreeNode::build(tripts,BoundingBox(bboxMin, bboxMax), 0, depth_max);

	/*
	 * Find the normals in case they are not present
	 */
	if(normals.size())
	{
		normals.resize(0);
		std::cout << "Normals are already defined! But still doing flat shading." << std::endl;
	}

	if(getpot.search("-s") && (triangles.size() > 200))
	{
		/*
		 * Smooth shading
		 */

		std::vector<double> vert_sigA(positions.size(), 0);
		std::vector<Vector> vert_sigNA(positions.size(), Vector(0.));

		// Find the sums of A and NA for each vertex
		for (auto it = triangles.begin() ; it != triangles.end() ; ++it)
		{
			Triangle &tri = *it;

			// Extract vertex positions from the mesh data.
			Vector const &p0 = positions[tri[0].pi];
			Vector const &p1 = positions[tri[1].pi];
			Vector const &p2 = positions[tri[2].pi];

			// Two of triangle edge
			const Vector e1 = p1 - p0;
			const Vector e2 = p2 - p0;

			// Find normal and area
			const Vector vecareatot = e1.cross(e2);
			const double Atot = vecareatot.length();

			for (int i = 0 ; i < 3  ; i++)
			{
				tri[i].ni = tri[i].pi;
				vert_sigA[tri[i].ni] += Atot;
				vert_sigNA[tri[i].ni] += vecareatot;
			}
		}

		// Find the smoothed normal at each vertex
		for (uint pi = 0 ; pi < positions.size() ; pi++)
		{
			normals.push_back( (1. / vert_sigA[pi]) * vert_sigNA[pi] );
		}
	} // End of smooth shading

	else
	{
		/*
		 * Flat shading
		 */
		for (auto it = triangles.begin() ; it != triangles.end() ; ++it)
		{
			Triangle &tri = *it;

			// Extract vertex positions from the mesh data.
			Vector const &p0 = positions[tri[0].pi];
			Vector const &p1 = positions[tri[1].pi];
			Vector const &p2 = positions[tri[2].pi];

			// Two of triangle edge
			const Vector e1 = p1 - p0;
			const Vector e2 = p2 - p0;

			// Find normal and area
			const Vector nsurf = e1.cross(e2).normalized();
			tri[0].ni = tri[1].ni = tri[2].ni = normals.size();
			normals.push_back(nsurf);

		}
	} // End of flat shading
}
