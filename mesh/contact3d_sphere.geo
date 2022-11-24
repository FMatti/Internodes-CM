l = 0.2;

// lower points
Point(1) = {0, 0, 0, l};
Point(2) = {1, 0, 0, l};
Point(3) = {1, 1, 0, l};
Point(4) = {0, 1, 0, l};
Point(11) = {0, 0, 1, l};
Point(12) = {1, 0, 1, l};
Point(13) = {1, 1, 1, l};
Point(14) = {0, 1, 1, l};

r = 0.5;
h = 1.4;

// circle center point
Point(21) = {0.5, h, 0.5, l};
Physical Point("blocked_nodes") = {21};

// upper points
Point(31) = {0.5, h, 0.5 - r, l};
Point(32) = {0.5, h, 0.5 + r, l};
Point(33) = {0.5 - r, h, 0.5, l};
Point(34) = {0.5 + r, h, 0.5, l};

Point(35) = {0.5, h - r, 0.5, l};

Circle(36) = {31, 21, 35};
Circle(37) = {32, 21, 35};
Circle(38) = {33, 21, 35};
Circle(39) = {34, 21, 35};
Circle(40) = {31, 21, 33};
Circle(41) = {33, 21, 32};
Circle(42) = {32, 21, 34};
Circle(43) = {34, 21, 31};


// lower lines
Line(1) = {1, 2};
Physical Curve("lower_bottom_front") = {1};
Line(2) = {2, 12};
Physical Curve("lower_bottom_left") = {2};
Line(3) = {12, 11};
Physical Curve("lower_bottom_back") = {3};
Line(4) = {11, 1};
Physical Curve("lower_bottom_right") = {4};

Line(5) = {1, 4};
Physical Curve("lower_front_right") = {5};
Line(6) = {2, 3};
Physical Curve("lower_front_left") = {6};
Line(7) = {12, 13};
Physical Curve("lower_back_left") = {7};
Line(8) = {11, 14};
Physical Curve("lower_back_right") = {8};

Line(9) = {3, 4};
Physical Curve("lower_top_front") = {9};
Line(10) = {4, 14};
Physical Curve("lower_top_right") = {10};
Line(11) = {14, 13};
Physical Curve("lower_top_back") = {11};
Line(12) = {13, 3};
Physical Curve("lower_top_left") = {12};

// lower surfaces
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("lower_bottom") = {1};
Curve Loop(2) = {-1, 5, -9, -6};
Plane Surface(2) = {2};
Physical Surface("lower_front") = {2};
Curve Loop(3) = {-2, 6, -12, -7};
Plane Surface(3) = {3};
Physical Surface("lower_left") = {3};
Curve Loop(4) = {-3, 7, -11, -8};
Plane Surface(4) = {4};
Physical Surface("lower_back") = {4};
Curve Loop(5) = {-4, 8, -10, -5};
Plane Surface(5) = {5};
Physical Surface("lower_right") = {5};
Curve Loop(6) = {9, 10, 11, 12};
Plane Surface(6) = {6};
Physical Surface("lower_top") = {6};


Curve Loop(11) = {-40, 36, -38};
Surface(11) = {11};
Curve Loop(12) = {-41, 38, -37};
Surface(12) = {12};
Curve Loop(13) = {-42, 37, -39};
Surface(13) = {13};
Curve Loop(14) = {-43, 39, -36};
Surface(14) = {14};
Physical Surface("upper_bottom") = {11, 12, 13, 14};
Curve Loop(15) = {40, 41, 42, 43};
Plane Surface(15) = {15};
Physical Surface("upper_top") = {15};

// lower volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Volume("body_lower") = {1};

// upper volume
Surface Loop(2) = {11, 12, 13, 14, 15};
Volume(2) = {2};
Physical Volume("body_upper") = {2};
