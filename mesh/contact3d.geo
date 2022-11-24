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

// upper points
Point(5) = {0, 1.1, 0, l};
Point(6) = {1, 1.1, 0, l};
Point(7) = {1, 2, 0, l};
Point(8) = {0, 2, 0, l};
Point(15) = {0, 1.1, 1, l};
Point(16) = {1, 1.1, 1, l};
Point(17) = {1, 2, 1, l};
Point(18) = {0, 2, 1, l};

// circle center points
Point(21) = {0.5, 2.1, 0, l};
Point(22) = {0.5, 2.1, 1, l};
Physical Point("blocked_nodes") = {21, 22};

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

// upper lines
Circle(21) = {5, 21, 6};
Physical Curve("upper_bottom_front") = {21};
Line(22) = {6, 16};
Physical Curve("upper_bottom_left") = {22};
Circle(23) = {16, 22, 15};
Physical Curve("upper_bottomback") = {23};
Line(24) = {15, 5};
Physical Curve("upper_bottom_right") = {24};

Line(25) = {5, 8};
Physical Curve("upper_front_right") = {25};
Line(26) = {6, 7};
Physical Curve("upper_front_left") = {26};
Line(27) = {16, 17};
Physical Curve("upper_back_left") = {27};
Line(28) = {15, 18};
Physical Curve("upper_back_right") = {28};

Line(29) = {7, 8};
Physical Curve("upper_top_front") = {29};
Line(30) = {8, 18};
Physical Curve("upper_top_right") = {30};
Line(31) = {18, 17};
Physical Curve("upper_top_back") = {31};
Line(32) = {17, 7};
Physical Curve("upper_top_left") = {32};

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

// upper surfaces
Curve Loop(11) = {24, 21, 22, 23};
Surface(11) = {11};
Physical Surface("upper_bottom") = {11};
Curve Loop(12) = {-21, 25, -29, -26};
Plane Surface(12) = {12};
Physical Surface("upper_front") = {12};
Curve Loop(13) = {-22, 26, -32, -27};
Plane Surface(13) = {13};
Physical Surface("upper_left") = {13};
Curve Loop(14) = {-23, 27, -31, -28};
Plane Surface(14) = {14};
Physical Surface("upper_back") = {14};
Curve Loop(15) = {-24, 28, -30, -25};
Plane Surface(15) = {15};
Physical Surface("upper_right") = {15};
Curve Loop(16) = {29, 30, 31, 32};
Plane Surface(16) = {16};
Physical Surface("upper_top") = {16};

// lower volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Volume("body_lower") = {1};

// upper surfaces
Surface Loop(2) = {11, 12, 13, 14, 15, 16};
Volume(2) = {2};
Physical Volume("body_upper") = {2};

// field
Field[1] = Box;
