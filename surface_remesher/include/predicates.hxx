/*
 * predicates.hxx
 *
 *  Created on: Dec 1, 2016
 *      Author: hooshi
 *
 *  Taken from TetGen, which in turn takes it from Jonathan Shewchuck.
 */

#ifndef INCLUDE_PREDICATES_HXX_
#define INCLUDE_PREDICATES_HXX_


#define REAL double

/*               Return a positive value if the points pa, pb, and pc occur  */
/*               in counterclockwise order; a negative value if they occur   */
/*               in clockwise order; and zero if they are collinear.  The    */
/*               result is also a rough approximation of twice the signed    */
/*               area of the triangle defined by the three points.           */
double orient2d(double* pa, double* pb, double* pc);
double orient2dfast(double* pa, double* pb, double* pc);
double orient2dslow(double* pa, double* pb, double* pc);


/*               Return a positive value if the point pd lies inside the     */
/*               circle passing through pa, pb, and pc; a negative value if  */
/*               it lies outside; and zero if the four points are cocircular.*/
/*               The points pa, pb, and pc must be in counterclockwise       */
/*               order, or the sign of the result will be reversed.          */
double incircle(double* pa, double* pb, double* pc, double* pd);
double incirclefast(double* pa, double* pb, double* pc, double* pd);
double incircleslow(double* pa, double* pb, double* pc, double* pd);

#endif /* INCLUDE_PREDICATES_HXX_ */
