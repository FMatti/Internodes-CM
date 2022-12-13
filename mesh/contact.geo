SetFactory("OpenCASCADE");
// Parameters
// Mesh parameter at point
lc=1.0;
// Radius of the sphere
R=0.5;
// Length
L=2;
// Height
H=2;



Point(1) = {-L, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {-L, H, 0, lc};
Point(5) = {0, H+R, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {-0, 2.5, 0, 0.5, 0, 2*Pi};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
Physical Surface("lower_body", 6) = {1};
//+
Physical Surface("upper_body", 7) = {2};
//+
Physical Curve("bottom_edge_lower_body", 8) = {1};
//+
Physical Curve("right_edge_lower_body", 9) = {2};
//+
Physical Curve("top_edge_lower_body", 10) = {3};
//+
Physical Curve("left_edge_lower_body", 11) = {4};
//+
Physical Curve("boundary_upper_body", 12) = {5};