// Inputs
sidex = 0.2;
sidey = 1;
sidez = 0.2;

cl = 0.05;

// Geometry
Point(1) = {0      , sidey     , 0    , cl};
Point(2) = {0      , sidey/2   , 0    , cl};
Point(3) = {0      , 0        , 0    , cl};
Point(4) = {sidex/2 , 0        , 0    , cl};
Point(5) = {sidex   , 0        , 0    , cl};
Point(6) = {sidex   , sidey/2   , 0    , cl};
Point(7) = {sidex   , sidey     , 0    , cl};
Point(8) = {sidex/2 , sidey     , 0    , cl};
 
Line(1) = {1, 2, 3};				// bottom line
Line(2) = {3, 4, 5};				// right line
Line(3) = {5, 6, 7};				// top line
Line(4) = {7, 8, 1};				// left line


Line Loop(1) = {1, 2, 3, 4}; 	
Plane Surface(1) = {1};

Transfinite Line {2,4} = 4;
Transfinite Line {1,3} = 21;
Transfinite Surface "*";

ext[] = Extrude {0, 0, -sidez}
{
Surface{1};
Layers{3};
};

Physical Surface(1) = {1, -ext[0],-ext[2],-ext[3],-ext[4],-ext[5]};
