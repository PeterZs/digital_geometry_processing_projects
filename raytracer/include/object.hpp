#ifndef OBJECT_H
#define OBJECT_H

#include "basic.hpp"
#include "GetPot"

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <memory>

// Read the command line
// THis is not very nice
extern GetPot getpot;

// The type of a "parameter list", e.g. mapping from strings to sets of numbers.
class TreeNode;
typedef std::map<std::string, std::vector<double> > ParamList;

struct BoundingBox
{
	Vector minp;
	Vector maxp;

	BoundingBox(){}
	BoundingBox(const Vector &min_in, const Vector &max_in):
		minp(min_in), maxp(max_in)
	{
		// Sanity check
		assert(minp[XX] < maxp[XX]);
		assert(minp[YY] < maxp[YY]);
		assert(minp[ZZ] < maxp[ZZ]);
	}

	// Break along longest axis
	void divide_longest_axis(BoundingBox &b1, BoundingBox &b2)
	{
		const double dx = maxp[XX] - minp[XX];
		const double dy = maxp[YY] - minp[YY];
		const double dz = maxp[ZZ] - minp[ZZ];

		b1.minp = minp;	b1.maxp = maxp;
		b2.minp = minp;	b2.maxp = maxp;

		if ( (dx > dy) && (dx > dz) )
		{
			b1.minp[XX] = b2.maxp[XX] = (minp[XX] + maxp[XX]) / 2.;
		}
		else if ( (dy > dx) && (dy > dz) )
		{
			b1.minp[YY] = b2.maxp[YY] = (minp[YY] + maxp[YY]) / 2.;
		}
		else
		{
			b1.minp[ZZ] = b2.maxp[ZZ] = (minp[ZZ] + maxp[ZZ]) / 2.;
		}
	}

	//Intersection with another bounding box
	// There are cases which this does not consider
	// But I think this cases should give reasonable performance.
	bool is_inside(const Vector &v) const
	{
		if ( (v[XX] > minp[XX]) && (v[XX] < maxp[XX]) &&
				(v[YY] > minp[YY]) && (v[YY] < maxp[YY])&&
				(v[ZZ] > minp[ZZ]) && (v[ZZ] < maxp[ZZ])	)
			return true;
		else
			return false;
	}

	bool is_containing(const BoundingBox &bb2) const
	{
		// something taken from the web
		return(
				(this->maxp[XX] > bb2.minp[XX]) &&
				(this->minp[XX]< bb2.maxp[XX]  )&&
				(this->maxp[YY] > bb2.minp[YY] )&&
				(this->minp[YY] < bb2.maxp[YY] )&&
				(this->maxp[ZZ] > bb2.minp[ZZ] )&&
				(this->minp[ZZ] < bb2.maxp[ZZ])
		);
	}
};

// A class to encapsulate all parameters for a material.
class Material
{
public:
	// Constructors
	Material() {};
	Material(ParamList &params) { init(params); }

	void init(ParamList &params)
	{
		#define SET_VECTOR(_name) _name = params[#_name];	
// !!!!! edited lines start
		SET_VECTOR(ambient)
		SET_VECTOR(diffuse)
		SET_VECTOR(specular)
// !!!!! edited lines end
		SET_VECTOR(emission)

		#define SET_FLOAT(_name) _name = params[#_name].size() ? params[#_name][0] : 0;
		SET_FLOAT(shininess)
		SET_FLOAT(shadow)
		SET_FLOAT(reflect)
	}

	// Ambient/diffuse/specular/emissive colors. 
	Vector ambient;
	Vector diffuse;
	Vector specular;	   
	Vector emission;

	// "Shininess" factor (specular exponent).
	double shininess;

	// Shadow coefficient, [0 -> no shadow, 1 -> black shadow]
	// everything in between is blended by that factor with the surface color
	double shadow;

	// Reflection coefficient [0 -> no reflection, 1 -> total reflection]
	// everything in between is blended by that factor with the surface color
	double reflect;
};


// Abstract base object class
class Object 
{
protected:
	static constexpr double inter_tol = 1e-14;
public:
    Matrix transform;   // Transformation from global to object space.
    Matrix i_transform; // Transformation from object to global space.
    Matrix n_transform; // Trasnformation to global space for normals.

    // Sets up the 3 transformations from the given global-to-object transform.
    void setup_transform(Matrix const &m)
	{
		transform = m;
		m.invert(i_transform);
		n_transform = i_transform.transpose();
	}

    // Intersect the object with the given ray in global space.
    // Returns true if there was an intersection, hit is updated with params.
    // NOTE: the function does not do anything with the initial value of hit.
    bool intersect(Ray ray, Intersection &hit) const;

    // Intersect the object with the given ray in object space.
    // This function is specific to each object subtype.
    // Returns true if there was an intersection, hit is updated with params.
    virtual bool localIntersect(Ray const &ray, Intersection &hit) const = 0;

    // Virtual Destructor
    // Mandatory in order to prevent memory leak.
    virtual ~Object(){}

    Material material; // This object's material.
};


// A sphere centred around the local origin with a certain radius.
class Sphere : public Object 
{
public:
    double radius;
    
    bool localIntersect(Ray const &ray, Intersection &hit) const;
};


// A plane at the origin using Z+ as the normal in object space.
class Plane : public Object 
{
public:
    virtual bool localIntersect(Ray const &ray, Intersection &hit) const;
};


// A conic about the Z+ axis, bounded along Z by zMin and zMax, 
// with radii radius1 and radius2.
class Conic : public Object 
{
public:
    double radius1, radius2;
    double zMin, zMax;

    bool localIntersect(Ray const &ray, Intersection &hit) const;
};


// A class to represent a single vertex of a polygon. The ints stored within
// are indices into the positions/texCoords/normals/colors vectors of the
// Object that it belongs to.
class Vertex {
public:
	// Indices into positions, texCoods, normals, and colors vectors.
	int pi, ti, ni, ci;

	Vertex() : pi(-1), ti(-1), ni(-1), ci(-1) {}

	Vertex(int pi, int ti, int ni, int ci) :
		pi(pi),
		ti(ti),
		ni(ni),
		ci(ci)
	{}
};


class Triangle
{
public:
	Vertex v[3];
	BoundingBox bbox;

	Triangle(Vertex const &v0, Vertex const &v1, Vertex const &v2)
	{
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
	}

	Vertex& operator[](int i) { return v[i]; }
	const Vertex& operator[](int i) const { return v[i]; }
};


class Mesh : public Object {
public:
	// Storage for positions/texCoords/normals/colors. Looked up by index.
	std::vector<Vector> positions;
	std::vector<Vector> texCoords;
	std::vector<Vector> normals;
	std::vector<Vector> colors;

	// Triangles are a triplet of vertices.
	std::vector<Triangle> triangles;

	// Bounding box
	Vector bboxMin, bboxMax;
	TreeNode* treeroot;

	// constructor
	Mesh(): treeroot(NULL){}
	~Mesh();

	// Read OBJ data from a given file.
	bool readOBJ(std::string const &filename);

	// Construct bounding box of vertex positions.
	// This must be called after mesh data has been initialized and before
	// raytracing begins.
	void updateBBox();

	// Intersections!
	bool localIntersect(Ray const &ray, Intersection &hit) const;
	bool localIntersect_bbox(Ray const &ray, TreeNode * subspace, Intersection &hit) const;

	// void init
	// setup kdtree for intersection tests, normals, ..
	void init();

private:
	// Compute the result of the implicit line equation in 2D 
	// for a given point and a line with the given endpoints.
	double implicitLineEquation(double p_x, double p_y,
		double e1_x, double e1_y,
		double e2_x, double e2_y) const;

	// Find the intersection point between the given ray and mesh triangle.
	// Return true iff an intersection exists, and fill the hit data
	// structure with information about it.
	bool intersectTriangle(Ray const &ray,
		Triangle const &tri,
		Intersection &hit) const;
};


#endif
