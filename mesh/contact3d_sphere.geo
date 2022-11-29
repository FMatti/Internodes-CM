l = 0.2; // Target mesh size 

lx = 2; // x-sidelength of cuboid
ly = 0.5; // y-sidelength of cuboid
lz = 2; // z-sidelength of cuboid

r = 0.5; // Radius of the circle
d = 0.05; // Initial penetration overlap of sphere

a = 0.35; // Maximum expected radius of contact area
n = 5; // Number of mesh refinements for candidates


// ---> PRIMARY

// Cuboid corners
Point(1) = {-lx/2, -ly, -lz/2, l};
Point(2) = {lx/2, -ly, -lz/2, l};
Point(3) = {lx/2, 0, -lz/2, l};
Point(4) = {-lx/2, 0, -lz/2, l};
Point(5) = {-lx/2, -ly, lz/2, l};
Point(6) = {lx/2, -ly, lz/2, l};
Point(7) = {lx/2, 0, lz/2, l};
Point(8) = {-lx/2, 0, lz/2, l};

Point(9) = {lx/2, 0, 0, l};
Point(10) = {0, 0, -lz/2, l};
Point(11) = {-lx/2, 0, 0, l};
Point(12) = {0, 0, lz/2, l};

// Cuboid edges
Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};
Line(5) = {1, 4};
Line(6) = {2, 3};
Line(7) = {6, 7};
Line(8) = {5, 8};
Line(9) = {3, 10};
Line(10) = {4, 11};
Line(11) = {8, 12};
Line(12) = {7, 9};

Line(13) = {10, 4};
Line(14) = {12, 7};
Line(15) = {9, 3};
Line(16) = {11, 8};

// Cuboid surfaces
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {-1, 5, -9, -13, -6};
Plane Surface(2) = {2};
Curve Loop(3) = {-2, 6, -15, -12, -7};
Plane Surface(3) = {3};
Curve Loop(4) = {-3, 7, -14, -11, -8};
Plane Surface(4) = {4};
Curve Loop(5) = {-4, 8, -16, -10, -5};
Plane Surface(5) = {5};

Point(13) = {0, 0, a, l};
Point(14) = {a, 0, 0, l};
Point(15) = {0, 0, -a, l};
Point(16) = {-a, 0, 0, l};

Point(17) = {0, 0, 0, l};

Circle(17) = {13, 17, 14};
Circle(18) = {14, 17, 15};
Circle(19) = {15, 17, 16};
Circle(20) = {16, 17, 13};
Line(21) = {9, 14};
Line(22) = {10, 15};
Line(23) = {11, 16};
Line(24) = {12, 13};

Delete { Point{17}; }

Curve Loop(6) = {15, 9, 22, -18, -21};
Plane Surface(6) = {6};

Curve Loop(7) = {13, 10, 23, -19, -22};
Plane Surface(7) = {7};

Curve Loop(8) = {16, 11, 24, -20, -23};
Plane Surface(8) = {8};

Curve Loop(9) = {14, 12, 21, -17, -24};
Plane Surface(9) = {9};

Curve Loop(10) = {17, 18, 19, 20};
Plane Surface(10) = {10};

Transfinite Curve{17} = n;
Transfinite Curve{18} = n;
Transfinite Curve{19} = n;
Transfinite Curve{20} = n;

Physical Curve("primary_transfinite") = {17, 18, 19, 20};

Physical Surface("primary_fixed") = {1};
Physical Surface("primary_bulk") = {2, 3, 4, 5, 6, 7, 8, 9};
Physical Surface("primary_candidates") = {10};

// Cuboid
S = news;
Surface Loop(S) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

V = newv;
Volume(V) = {S};
Physical Volume("primary") = {V};


// ---> SECONDARY

Geometry.ExtrudeSplinePoints = 10;

// Opening angle for candidate nodes
phi = Asin(a/r); 

// Sphere center point
C = newp;
Point(C) = {0, r-d, 0, l};

// Sphere bottom point
P0 = newp;
Point(P0) = {0, -d, 0, l};

// Sphere candidates end-point
P1 = newp;
Point(P1) = {r*Sin(phi), -r*Cos(phi) + r-d, 0, l};

// Sphere end-point
P2 = newp;
Point(P2) = {r, r-d, 0, l};

// Circle connecting bottom point to candidates end-point
L1 = newl;
Circle(L1) = {P0, C, P1};

// Circle connecting candidates end-point to end-point
L2 = newl;
Circle(L2) = {P1, C, P2};

// Line connecting end-point to sphere's center point
L3 = newl;
Line(L3) = {P2, C};

// Extrude candidates surface
E1a[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{L1}; };
E1b[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{E1a[0]}; };

// Extrude bulk surface
E2a[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{L2}; };
E2b[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{E2a[0]}; };

// Extrude top surface
E3a[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{L3}; };
E3b[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{E3a[0]}; };

Transfinite Curve(27) = n;
Transfinite Curve(30) = n;
Transfinite Curve(32) = 2*n;
Transfinite Curve(36) = 2*n;

Physical Curve("secondary_transfinite") = {27, 30, 32, 36};

Physical Surface("secondary_candidates") = {E1a[1], E1b[1]};
Physical Surface("secondary_bulk") = {E2a[1], E2b[1]};
Physical Surface("secondary_fixed") = {E3a[1], E3b[1]};

// Semi-sphere
S = news;
Surface Loop(S) = {E1a[1], E1b[1], E2a[1], E2b[1], E3a[1], E3b[1]};

V = newv;
Volume(V) = {S};
Physical Volume("secondary") = {V};
