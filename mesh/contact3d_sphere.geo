l = 0.1; // Target mesh size 

lx = 1; // x-sidelength of cuboid
ly = 0.5; // y-sidelength of cuboid
lz = 1; // z-sidelength of cuboid

r = 0.5; // Radius of the sphere
d = 0.05; // Initial penetration overlap of sphere
phi = Pi/4; // Opening angle for candidate nodes


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

// Cuboid edges
Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};
Line(5) = {1, 4};
Line(6) = {2, 3};
Line(7) = {6, 7};
Line(8) = {5, 8};
Line(9) = {3, 4};
Line(10) = {4, 8};
Line(11) = {8, 7};
Line(12) = {7, 3};

// Cuboid surfaces
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {-1, 5, -9, -6};
Plane Surface(2) = {2};
Curve Loop(3) = {-2, 6, -12, -7};
Plane Surface(3) = {3};
Curve Loop(4) = {-3, 7, -11, -8};
Plane Surface(4) = {4};
Curve Loop(5) = {-4, 8, -10, -5};
Plane Surface(5) = {5};
Curve Loop(6) = {9, 10, 11, 12};
Plane Surface(6) = {6};

Physical Surface("primary_fixed") = {1};
Physical Surface("primary_bulk") = {2, 3, 4, 5};
Physical Surface("primary_candidates") = {6};

// Cuboid
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Volume("primary") = {1};


// ---> SECONDARY

Geometry.ExtrudeSplinePoints = 10;

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

Physical Surface("secondary_candidates") = {E1a[1], E1b[1]};
Physical Surface("secondary_bulk") = {E2a[1], E2b[1]};
Physical Surface("secondary_fixed") = {E3a[1], E3b[1]};

// Semi-sphere
S = news;
Surface Loop(S) = {E1a[1], E1b[1], E2a[1], E2b[1], E3a[1], E3b[1]};

V = newv;
Volume(V) = {S};
Physical Volume("secondary") = {V};
