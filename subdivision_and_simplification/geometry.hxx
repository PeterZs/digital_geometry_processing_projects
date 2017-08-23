/*
 * edit_mesh.hxx
 *
 * Header file for the class EditMesh, which encapsulates a manifold
 * triangular mesh.
 */

#ifndef GEOMETRY_HXX
#define GEOMETRY_HXX

#include <cmath>
#include <iostream>

#include <Eigen/Geometry>


#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define ABS(x) ((x) > 0 ? (x) : (-(x)) )

typedef std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > VecOfVerts;

namespace geo
{

    using Eigen::Vector3d;

    struct Plane
    {
        Vector3d b;
        Vector3d n;
    };

    inline double dist_from_line(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3)
    {
        const Vector3d e23 = v3 - v2;
        const Vector3d n23 = e23 / sqrt(e23.dot(e23));
        const Vector3d e21 = v1 - v2;
        const Vector3d e2n = e21 - (e21.dot(n23) * n23); 
        const double ans = sqrt(e2n.dot(e2n));
        //if(!finite(ans)) std::cout << e21.dot(e21) - e21.dot(n23) << std::endl;
        assert(finite(ans));
        return ans;
    }
    
    inline Vector3d normal(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3)
    {
        Vector3d e1 = v2 - v1;
        Vector3d e2 = v3 - v1;
        Vector3d ans = e1.cross(e2);
        // std::cout << "v1: " << v1.transpose() << std::endl;
        // std::cout << "v2: " << v2.transpose() << std::endl;
        // std::cout << "v3: " << v3.transpose() << std::endl;
                
        // std::cout << "e1: " << e1.transpose() << std::endl;
        // std::cout << "e2: " << e2.transpose() << std::endl;
        // std::cout << "normal: " << ans.transpose() << std::endl;
        // assert(finite(ans.x()));
        // assert(finite(ans.y()));
        // assert(finite(ans.z()));            
        
        return ans / sqrt( ans.dot(ans) );
    }

 
    inline double angle_n(const Vector3d& n1, const Vector3d& n2)
    {
        const double cs = n1.dot(n2);
        const double tol = 1e-3;
        double ans;
        if((cs > 1) && (cs < 1+tol)) ans = 0;
        else if ( (cs < -1) && (cs > -1-tol) ) ans = -M_PI;
        else  ans = acos(cs);
        // std::cout << n1.transpose() << std::endl;
        // std::cout << n2.transpose() << std::endl;
        // std::cout << n1.dot(n2) << std::endl;
        
        // printf("%lf \n", ans);
        if(!finite(ans))
        {
            printf("bade n1.dot(n2) %.8lf \n" , n1.dot(n2) );
        }
        assert(finite(ans));
        return ans;
    }

    inline void avg_plane(VecOfVerts& xyz, const std::size_t v, const std::vector<std::size_t>& rv, Plane& ans)
    {
        Vector3d sumcent = Vector3d::Zero(), sumn = Vector3d::Zero(), cent, normal, e1, e2, cross;
        double area, sumarea = 0;

        // for (uint i = 0 ; i < rv.size() ; i++)
        // {
        //     std::cout << rv[i] << "-> " << xyz[rv[i]].transpose() << std::endl;
        // }        
        // std::cout << std::endl;
        
        for (uint i = 0 ; i < rv.size() ; i++)
        {
            const std::size_t a = rv[i];
            const std::size_t b = rv[(i+1) % rv.size()];

            e1 = xyz[a] - xyz[v];
            e2 = xyz[b] - xyz[v];

            cross = e1.cross(e2);
            area = sqrt( cross.dot(cross) ) / 2.;
            normal = cross / (2. * area);
            cent = 1/3. * (xyz[a] + xyz[b] + xyz[v]);

            sumn += normal * area;
            sumcent += cent * area;
            sumarea += area;

            //debug:
            // std::cout << "xyz[0]: " << xyz[0].transpose() << std::endl;
            // std::cout << "xyz[a]: " << xyz[a].transpose() << std::endl;
            // std::cout << "xyz[v]: " << xyz[v].transpose() << std::endl; 
            // std::cout << "xyz[b]: " << xyz[b].transpose() << std::endl; 
     
            // std::cout << "area: " << area << " ";
            // std::cout << "center: " << cent.transpose() << std::endl;
        }

        ans.b = sumcent / sumarea;
        ans.n = sumn / sumarea;
        const double nsize = sqrt(ans.n.dot(ans.n));
        ans.n = ans.n / nsize;
        
        //assert( std::abs(nsize-1) < 1e-6 );
        // std::cout << "nsize ave normal: " << nsize << std::endl;
    }

    inline double dist_from_plane(const Plane& plane, const Vector3d& p)
    {
        return (p - plane.b).dot(plane.n);
    }

    inline void project_on_plane(const Plane& plane, const Vector3d& p, Vector3d& proj)
    {
        const double dist = dist_from_plane(plane, p);
        proj =  (p - dist * plane.n);
    }

    inline void div_plane(const Plane& aveplane, const Vector3d& v1, const Vector3d& v2, Plane& ans)
    {
        const Vector3d e = v2 - v1;
        Vector3d cross = e.cross(aveplane.n);
        ans.n = cross / sqrt( cross.dot(cross) );
        ans.b = v1;
    }

    // Invert a 3x3 Matrix
    // Written by C. Ollivier-Gooch
    inline bool Invert3x3(const double Block[3][3], double Inverse[3][3], const double tol=1e-6)
    {
        const double Det = (+ Block[0][0]*Block[1][1]*Block[2][2]
                            + Block[0][1]*Block[1][2]*Block[2][0]
                            + Block[0][2]*Block[1][0]*Block[2][1]
                            - Block[0][2]*Block[1][1]*Block[2][0]
                            - Block[0][1]*Block[1][0]*Block[2][2]
                            - Block[0][0]*Block[1][2]*Block[2][1]);
        if (std::abs(Det) < tol) return false;
        const double DetInv = 1/Det;
        
        // Expand by minors to compute the inverse.
        Inverse[0][0] = + DetInv * (Block[1][1]*Block[2][2] -
                                    Block[2][1]*Block[1][2]); 
        Inverse[1][0] = - DetInv * (Block[1][0]*Block[2][2] -
                                    Block[2][0]*Block[1][2]); 
        Inverse[2][0] = + DetInv * (Block[1][0]*Block[2][1] -
                                    Block[2][0]*Block[1][1]); 
        Inverse[0][1] = - DetInv * (Block[0][1]*Block[2][2] -
                                    Block[2][1]*Block[0][2]); 
        Inverse[1][1] = + DetInv * (Block[0][0]*Block[2][2] -
                                    Block[2][0]*Block[0][2]); 
        Inverse[2][1] = - DetInv * (Block[0][0]*Block[2][1] -
                                    Block[2][0]*Block[0][1]); 
        Inverse[0][2] = + DetInv * (Block[0][1]*Block[1][2] -
                                    Block[1][1]*Block[0][2]); 
        Inverse[1][2] = - DetInv * (Block[0][0]*Block[1][2] -
                                    Block[1][0]*Block[0][2]); 
        Inverse[2][2] = + DetInv * (Block[0][0]*Block[1][1] -
                                    Block[1][0]*Block[0][1]);
        return true;
    }

    inline void MultVec(double A[3][3], const double Vec[3], double Result[3])
    {
        Result[0] = A[0][0]*Vec[0] + A[0][1]*Vec[1] + A[0][2]*Vec[2]; 
        Result[1] = A[1][0]*Vec[0] + A[1][1]*Vec[1] + A[1][2]*Vec[2]; 
        Result[2] = A[2][0]*Vec[0] + A[2][1]*Vec[1] + A[2][2]*Vec[2]; 
    }

}
#endif /*GEOMETRY_HXX*/
