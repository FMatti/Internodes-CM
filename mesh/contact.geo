l = 0.2;
//+
Point(1) = {0, 0, 0, l};
Point(2) = {1, 0, 0, l};
Point(3) = {1, 1, 0, l};
Point(4) = {0, 1, 0, l};
Point(5) = {0, 1.1, 0, l};
Point(6) = {1, 1.1, 0, l};
Point(7) = {1, 2, 0, l};
Point(8) = {0, 2, 0, l};
Point(9) = {0.5, 2.1, 0, l};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {5, 9, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Physical Surface("body_lower") = {1};
//+
Physical Surface("body_upper") = {2};
//+
Physical Curve("lower_bottom") = {1};
//+
Physical Curve("lower_right") = {2};
//+
Physical Curve("lower_left") = {4};
//+
Physical Curve("lower_top") = {3};
//+
Physical Curve("upper_bottom") = {5};
//+
Physical Curve("upper_right") = {6};
//+
Physical Curve("upper_left") = {8};
//+
Physical Curve("upper_top") = {7};
//+
Physical Point("blocked_nodes") = {1};
//+
Field[1] = Box;
