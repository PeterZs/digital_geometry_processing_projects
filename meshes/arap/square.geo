// Inputs
side = 1;
cl = 0.05;

// Geometry
Point(1) = {0      , side     , 0    , cl};
Point(2) = {0      , side/2   , 0    , cl};
Point(3) = {0      , 0        , 0    , cl};
Point(4) = {side/2 , 0        , 0    , cl};
Point(5) = {side   , 0        , 0    , cl};
Point(6) = {side   , side/2   , 0    , cl};
Point(7) = {side   , side     , 0    , cl};
Point(8) = {side/2 , side     , 0    , cl};
Point(9) = {side/2 , side/2   , 0    , cl/20};

 
Line(1) = {1, 2, 3};				// bottom line
Line(2) = {3, 4, 5};				// right line
Line(3) = {5, 6, 7};				// top line
Line(4) = {7, 8, 1};				// left line


Line Loop(1) = {1, 2, 3, 4}; 	
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Point {9} In Surface {1};
