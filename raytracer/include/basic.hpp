#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cassert>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef ABS_FLOAT
#define ABS_FLOAT(X) ((X)>0?(X):(-(X)))
#endif

inline double rad2deg(double rad)
{
	return rad * 360.0 / (2 * PI);
}

inline double deg2rad(double deg)
{
	return deg * 2 * PI / 360.0;
}

// Enums
enum Direction{XX=0, YY, ZZ};

// A class to represent and perform basic math operations on homogeneous
// coordinates, RGB colors, etc..
class Vector 
{
private:
    double data[4];

public:

    // Direct and copy constructors.
    Vector(double x=0.0, double y=0.0, double z=0.0, double w=0.0)
	{
		data[0] = x;
		data[1] = y;
		data[2] = z;
		data[3] = w;
	}
	Vector(std::vector<double> const &vec)
	{
		for (size_t i = 0; i < 4; i++) {
			data[i] = vec.size() > i ? vec[i] : 0;
		}
	}
	Vector(Vector const &other)
	{
		for (int i = 0; i < 4; i++) {
			data[i] = other.data[i];
		}
	}
    
    // Prints a vector like "(x,y,z,w)", without a newline.
	void print() const
	{
		std::cout << '(';
		for (int i = 0; i < 4; i++) {
			std::cout << (i ? "," : "") << data[i];
		}
		std::cout << ')';
	}
    
    // Indexing operators. Allows for getting and setting components via
    // array syntax.
    const double& operator[](int index) const { return data[index]; }
    double& operator[](int index) { return data[index]; }
    
    // Scalar multiplication of all components.
    // E.g.: `Vector b = a * k` where a is a Vector, and k is a double.
	Vector operator*(double value) const
	{
		Vector result;
		for (int i = 0; i < 4; i++) {
			result[i] = data[i] * value;
		}
		return result;
	}
	friend Vector operator*(double scale, const Vector& vec)
	{
		return vec * scale;
	}
    
    // Scalar division of all components.
    // E.g.: `Vector b = a / k` where a is a Vector, and k is a double.
	Vector operator/(double value) const
	{
		Vector result;
		for (int i = 0; i < 4; i++) {
			result[i] = data[i] / value;
		}
		return result;
	}

    // Component-wise multiplation of two vectors.
    // E.g.: `Vector c = a * b` performs `c[i] = a[i] * b[i]`.
	Vector operator*(Vector const &other) const
	{
		Vector result;
		for (int i = 0; i < 4; i++) {
			result[i] = data[i] * other.data[i];
		}
		return result;
	}

	const Vector operator-(void) const
	{
		return Vector(-data[0], -data[1], -data[2], data[3]);
	}

	bool operator!=(const Vector& rhs) const
	{
		if ((*this - rhs).length2() < 1e-12)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	bool operator==(const Vector& rhs) const
	{
		if (*this != rhs)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
    
    // ======================================================================
    // Every method below this line ignores the 4th coordinate. Ergo, if your
    // W is not 1, then these will return incorrect homogeneous results.
    // ======================================================================

    // Length and squared length of the represented coordinate.
	double length() const
	{
		double sum = 0;
		// Ignoring W.
		for (int i = 0; i < 3; i++) {
			sum += data[i] * data[i];
		}
		return sqrt(sum);
	}
	double length2() const
	{
		double sum = 0;
		// Ignoring W.
		for (int i = 0; i < 3; i++) {
			sum += data[i] * data[i];
		}
		return sum;
	}

    // Normalize represented coordinate (scale to length 1) in place; modifies the vector.
	void normalize()
	{
		double len = length();
		// Ignoring W.

		if (!(ABS_FLOAT(len - 1.0) < 1e-6))
		{
			for (int i = 0; i < 3; i++) {
				data[i] /= len;
			}
		}
	}

    // Get a normalized copy of the represented coordinate.
	Vector normalized() const
	{
		Vector v(*this);
		v.normalize();
		return v;
	}
    
    // Addition and subtraction of represented coordinates.
	Vector operator+(Vector const &other) const
	{
		Vector result;
		// Ignoring W.
		for (int i = 0; i < 3; i++) {
			result[i] = data[i] + other.data[i];
		}
		return result;
	}
	Vector operator-(Vector const &other) const
	{
		Vector result;
		// Ignoring W.
		for (int i = 0; i < 3; i++) {
			result[i] = data[i] - other.data[i];
		}
		return result;
	}

	const Vector& operator+=(const Vector& rhs)
	{
		*this = *this + rhs;
		return *this;
	}
    
    // Dot product of two vectors. Optional flag indicates if we should
    // be including the 4th coordinate.
    // E.g.: `Vector c = a.dot(b)` for `c = a . b`.
	double dot(Vector const &other, bool homogeneous = false) const
	{
		double result = 0;
		// Ignoring W.
		for (int i = 0; i < (homogeneous ? 4 : 3); i++) {
			result += data[i] * other.data[i];
		}
		return result;
	}

    // Cross product of two vectors.
    // E.g.: `Vector c = a.cross(b)` for `c = a x b`. 
	Vector cross(Vector const &other, double w = 1.0) const
	{
		return Vector(
			data[1] * other.data[2] - other.data[1] * data[2],
			-data[0] * other.data[2] + other.data[0] * data[2],
			data[0] * other.data[1] - other.data[0] * data[1],
			w
			);
	}
};


inline std::ostream& operator<<(std::ostream &out, Vector const &vec)
{
	out << vec[0] << " " << vec[1] << " " << vec[2] << " " <<vec[3];
	return out;
}

// A class to represent 4x4 transformation matrices. They are stored in
// row-major order, so they are indexed by row and then column.
class Matrix 
{
private:
    Vector data[4];

public:    
    // Copy constructor.
	Matrix(Matrix const &other)
	{
		for (int i = 0; i < 4; i++) {
			data[i] = other.data[i];
		}
	}

    // Constructo from rows.
	Matrix(Vector const &a, Vector const &b, Vector const &c, Vector const &d)
	{
		data[0] = a;
		data[1] = b;
		data[2] = c;
		data[3] = d;
	}

    // Construct from 16 doubles.
	Matrix(double a = 0, double b = 0, double c = 0, double d = 0,
		double e = 0, double f = 0, double g = 0, double h = 0,
		double i = 0, double j = 0, double k = 0, double l = 0,
		double m = 0, double n = 0, double o = 0, double p = 1)
	{
		data[0] = Vector(a, b, c, d);
		data[1] = Vector(e, f, g, h);
		data[2] = Vector(i, j, k, l);
		data[3] = Vector(m, n, o, p);
	}
    
    // Build an identity matrix.
	static Matrix identity()
	{
		return Matrix(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
			);
	}

    // Build a matrix for rotating about an arbitrary axis.
	static Matrix rotation(double angle, Vector const &axis)
	{
		double c_angle = cos(angle);
		double v_angle = 1.0 - c_angle;
		double s_angle = sin(angle);

		return Matrix(
			axis[0] * axis[0] * v_angle + c_angle,
			axis[1] * axis[0] * v_angle - axis[2] * s_angle,
			axis[2] * axis[0] * v_angle + axis[1] * s_angle,
			0,
			axis[0] * axis[1] * v_angle + axis[2] * s_angle,
			axis[1] * axis[1] * v_angle + c_angle,
			axis[2] * axis[1] * v_angle - axis[0] * s_angle,
			0,
			axis[0] * axis[2] * v_angle - axis[1] * s_angle,
			axis[1] * axis[2] * v_angle + axis[0] * s_angle,
			axis[2] * axis[2] * v_angle + c_angle,
			0,
			0, 0, 0, 1
			);
	}
    
    // Build a matrix for translating.
	static Matrix translation(double x, double y, double z)
	{
		return Matrix(1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1);
	}

    // Build a matrix for scaling.
	static Matrix scale(double x, double y, double z)
	{
		return Matrix(x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1);
	}

    // Print out a matrix (split onto 4 lines if "pretty").
	void print(bool pretty = false) const
	{
		std::cout << '(';
		for (int r = 0; r < 4; r++) {
			std::cout << (r ? ";" : "") << (r && pretty ? "\n " : "");
			for (int c = 0; c < 4; c++) {
				std::cout << (c ? "," : "") << data[r][c];
			}
		}
		std::cout << ')';
	}
    
    // Indexing operators. Allows for getting/setting rows of a matrix,
    // which can also be indexed themselves.
    // E.g.: `m[0]` is the first row.
    //       `m[1][2]` is the 3rd element of the 2nd row.
    inline Vector  operator[](int index) const { return data[index]; }
    inline Vector& operator[](int index)       { return data[index]; }
    
    // Return a transposed copy of the matrix.
	Matrix transpose() const
	{
		Matrix result;
		for (int r = 0; r < 4; r++) {
			for (int c = 0; c < 4; c++) {
				result[r][c] = data[c][r];
			}
		}
		return result;
	}

    // Invert a matrix.
	bool invert(Matrix &inv) const
	{
		inv[0][0] = data[1][1] * data[2][2] * data[3][3] -
			data[1][1] * data[2][3] * data[3][2] -
			data[2][1] * data[1][2] * data[3][3] +
			data[2][1] * data[1][3] * data[3][2] +
			data[3][1] * data[1][2] * data[2][3] -
			data[3][1] * data[1][3] * data[2][2];
		inv[1][0] = -data[1][0] * data[2][2] * data[3][3] +
			data[1][0] * data[2][3] * data[3][2] +
			data[2][0] * data[1][2] * data[3][3] -
			data[2][0] * data[1][3] * data[3][2] -
			data[3][0] * data[1][2] * data[2][3] +
			data[3][0] * data[1][3] * data[2][2];
		inv[2][0] = data[1][0] * data[2][1] * data[3][3] -
			data[1][0] * data[2][3] * data[3][1] -
			data[2][0] * data[1][1] * data[3][3] +
			data[2][0] * data[1][3] * data[3][1] +
			data[3][0] * data[1][1] * data[2][3] -
			data[3][0] * data[1][3] * data[2][1];
		inv[3][0] = -data[1][0] * data[2][1] * data[3][2] +
			data[1][0] * data[2][2] * data[3][1] +
			data[2][0] * data[1][1] * data[3][2] -
			data[2][0] * data[1][2] * data[3][1] -
			data[3][0] * data[1][1] * data[2][2] +
			data[3][0] * data[1][2] * data[2][1];
		inv[0][1] = -data[0][1] * data[2][2] * data[3][3] +
			data[0][1] * data[2][3] * data[3][2] +
			data[2][1] * data[0][2] * data[3][3] -
			data[2][1] * data[0][3] * data[3][2] -
			data[3][1] * data[0][2] * data[2][3] +
			data[3][1] * data[0][3] * data[2][2];
		inv[1][1] = data[0][0] * data[2][2] * data[3][3] -
			data[0][0] * data[2][3] * data[3][2] -
			data[2][0] * data[0][2] * data[3][3] +
			data[2][0] * data[0][3] * data[3][2] +
			data[3][0] * data[0][2] * data[2][3] -
			data[3][0] * data[0][3] * data[2][2];
		inv[2][1] = -data[0][0] * data[2][1] * data[3][3] +
			data[0][0] * data[2][3] * data[3][1] +
			data[2][0] * data[0][1] * data[3][3] -
			data[2][0] * data[0][3] * data[3][1] -
			data[3][0] * data[0][1] * data[2][3] +
			data[3][0] * data[0][3] * data[2][1];
		inv[3][1] = data[0][0] * data[2][1] * data[3][2] -
			data[0][0] * data[2][2] * data[3][1] -
			data[2][0] * data[0][1] * data[3][2] +
			data[2][0] * data[0][2] * data[3][1] +
			data[3][0] * data[0][1] * data[2][2] -
			data[3][0] * data[0][2] * data[2][1];
		inv[0][2] = data[0][1] * data[1][2] * data[3][3] -
			data[0][1] * data[1][3] * data[3][2] -
			data[1][1] * data[0][2] * data[3][3] +
			data[1][1] * data[0][3] * data[3][2] +
			data[3][1] * data[0][2] * data[1][3] -
			data[3][1] * data[0][3] * data[1][2];
		inv[1][2] = -data[0][0] * data[1][2] * data[3][3] +
			data[0][0] * data[1][3] * data[3][2] +
			data[1][0] * data[0][2] * data[3][3] -
			data[1][0] * data[0][3] * data[3][2] -
			data[3][0] * data[0][2] * data[1][3] +
			data[3][0] * data[0][3] * data[1][2];
		inv[2][2] = data[0][0] * data[1][1] * data[3][3] -
			data[0][0] * data[1][3] * data[3][1] -
			data[1][0] * data[0][1] * data[3][3] +
			data[1][0] * data[0][3] * data[3][1] +
			data[3][0] * data[0][1] * data[1][3] -
			data[3][0] * data[0][3] * data[1][1];
		inv[3][2] = -data[0][0] * data[1][1] * data[3][2] +
			data[0][0] * data[1][2] * data[3][1] +
			data[1][0] * data[0][1] * data[3][2] -
			data[1][0] * data[0][2] * data[3][1] -
			data[3][0] * data[0][1] * data[1][2] +
			data[3][0] * data[0][2] * data[1][1];
		inv[0][3] = -data[0][1] * data[1][2] * data[2][3] +
			data[0][1] * data[1][3] * data[2][2] +
			data[1][1] * data[0][2] * data[2][3] -
			data[1][1] * data[0][3] * data[2][2] -
			data[2][1] * data[0][2] * data[1][3] +
			data[2][1] * data[0][3] * data[1][2];
		inv[1][3] = data[0][0] * data[1][2] * data[2][3] -
			data[0][0] * data[1][3] * data[2][2] -
			data[1][0] * data[0][2] * data[2][3] +
			data[1][0] * data[0][3] * data[2][2] +
			data[2][0] * data[0][2] * data[1][3] -
			data[2][0] * data[0][3] * data[1][2];
		inv[2][3] = -data[0][0] * data[1][1] * data[2][3] +
			data[0][0] * data[1][3] * data[2][1] +
			data[1][0] * data[0][1] * data[2][3] -
			data[1][0] * data[0][3] * data[2][1] -
			data[2][0] * data[0][1] * data[1][3] +
			data[2][0] * data[0][3] * data[1][1];
		inv[3][3] = data[0][0] * data[1][1] * data[2][2] -
			data[0][0] * data[1][2] * data[2][1] -
			data[1][0] * data[0][1] * data[2][2] +
			data[1][0] * data[0][2] * data[2][1] +
			data[2][0] * data[0][1] * data[1][2] -
			data[2][0] * data[0][2] * data[1][1];

		double det = data[0][0] * inv[0][0] +
			data[0][1] * inv[1][0] +
			data[0][2] * inv[2][0] +
			data[0][3] * inv[3][0];

		if (det == 0) {
			return false;
		}

		det = 1.0 / det;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				inv[i][j] *= det;
			}
		}

		return true;
	}

	//calculate leftUp 3x3 matrix's inverse matrix, useful for solving linear system
	bool get33Invert(Matrix &inv) const
	{
		const Matrix A = this->transpose();
		double determinant = +A(0, 0)*(A(1, 1)*A(2, 2) - A(2, 1)*A(1, 2))
			- A(0, 1)*(A(1, 0)*A(2, 2) - A(1, 2)*A(2, 0))
			+ A(0, 2)*(A(1, 0)*A(2, 1) - A(1, 1)*A(2, 0));
		if (ABS_FLOAT(determinant) < 1e-6)
		{
			return false;
		}
		else
		{
			Matrix &result = inv;
			double invdet = 1 / determinant;
			result(0, 0) = (A(1, 1)*A(2, 2) - A(2, 1)*A(1, 2))*invdet;
			result(1, 0) = -(A(0, 1)*A(2, 2) - A(0, 2)*A(2, 1))*invdet;
			result(2, 0) = (A(0, 1)*A(1, 2) - A(0, 2)*A(1, 1))*invdet;
			result(0, 1) = -(A(1, 0)*A(2, 2) - A(1, 2)*A(2, 0))*invdet;
			result(1, 1) = (A(0, 0)*A(2, 2) - A(0, 2)*A(2, 0))*invdet;
			result(2, 1) = -(A(0, 0)*A(1, 2) - A(1, 0)*A(0, 2))*invdet;
			result(0, 2) = (A(1, 0)*A(2, 1) - A(2, 0)*A(1, 1))*invdet;
			result(1, 2) = -(A(0, 0)*A(2, 1) - A(2, 0)*A(0, 1))*invdet;
			result(2, 2) = (A(0, 0)*A(1, 1) - A(1, 0)*A(0, 1))*invdet;

			return true;
		}
	}

	const double& operator()(int i, int j) const
	{
		return this->data[i][j];
	}
	double& operator()(int i, int j)
	{
		return this->data[i][j];
	}

    // Multiplication of a column Vector by the matrix.
    // E.g.: `Vector b = m * a` where m is a Matrix and a is a Vector.
	Vector operator*(Vector const &other) const
	{
		return Vector(
			other.dot(data[0], true),
			other.dot(data[1], true),
			other.dot(data[2], true),
			other.dot(data[3], true)
			);
	}
	Vector mult33(Vector const &other) const
	{
		return Vector(
			other.dot(data[0], false),
			other.dot(data[1], false),
			other.dot(data[2], false),
			other.dot(data[3], false)
			);
	}

    // Matrix multiplication.
    // E.g.: `Matrix c = a * b`.
	Matrix operator*(Matrix const &other) const
	{
		Matrix transposed = other.transpose();
		Matrix result;
		for (int r = 0; r < 4; r++) {
			for (int c = 0; c < 4; c++) {
				result[r][c] = data[r].dot(transposed.data[c], true);
			}
		}
		return result;
	}

    // In-place matrix multiplication.
    // E.g.: `a *= b` is the same as `a = a * b`.
	void operator*=(Matrix const &other)
	{
		Matrix copy = (*this) * other;
		for (int i = 0; i < 4; i++) {
			data[i] = copy.data[i];
		}
	}
};


// A class to encapsulate all of the information relating to a ray intersection.
class Intersection {
public:
	// How far along the ray the intersection occurred.
	double depth;

	// Location of the intersection.
	Vector position;

	// Surface normal at the point of intersection.
	Vector normal;

	Intersection(void) : depth(DBL_MAX), position(0), normal(0) {}
	Intersection(const double depth_in) : depth(depth_in), position(0), normal(0) {}
};


// A class to encapsulate a ray.
class Ray 
{
public:
	Ray(void) : origin(0, 0, 0), direction(0, 0, 1) {}
	Ray(double x, double y, double z, double dx, double dy, double dz) :
		origin(x, y, z), direction(dx, dy, dz) 
	{
	}

	Ray(Vector const &origin_, Vector const &direction_) :
		origin(origin_), direction(direction_)
	{
		origin[3] = 1; direction[3] = 0;
	}

	Vector origin;    // Ray origin
	Vector direction; // Ray direction, its length matters
	//if an intersection happens at origin + t * direction, then t is the depth
};


#endif
