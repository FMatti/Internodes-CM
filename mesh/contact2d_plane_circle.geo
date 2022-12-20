l = 0.3; // Target mesh size

lx = 2; // x-sidelength of cuboid
ly = 0.5; // y-sidelength of cuboid

r = 0.5; // Radius of the circle
d = 0.05; // Initial penetration overlap of sphere

a = 0.4; // Maximum expected radius of contact area
n = 20; // Number of mesh refinements for candidates


// ---> PRIMARY

// Nodes
Point(1) = {-lx/2, -ly, 0, l};
Point(2) = {lx/2, -ly, 0, l};
Point(3) = {lx/2, 0, 0, l};
Point(4) = {a, 0, 0, l};
Point(5) = {-a, 0, 0, l};
Point(6) = {-lx/2, 0, 0, l};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Physical Curve("primary_fixed") = {1};
Physical Curve("primary_bulk") = {2, 3, 5, 6};
Physical Curve("primary_candidates") = {4};

Transfinite Curve{4} = n;

// Surface
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
Physical Surface("primary") = {1};


// ---> SECONDARY

// Opening angle for candidate nodes
phi = Asin(a/r);

// Semi-circle center point
C = newp;
Point(C) = {0, r-d, 0, l};

// Semi-circle left corner point
P1 = newp;
Point(P1) = {-r, r-d, 0, l};

// Semi-circle candidate left end-point
P2 = newp;
Point(P2) = {-r*Sin(phi), -r*Cos(phi) + r-d, 0, l};

// Circle candidate right end-point
P3 = newp;
Point(P3) = {r*Sin(phi), -r*Cos(phi) + r-d, 0, l};

// Semi-circle right corner point
P4 = newp;
Point(P4) = {r, r-d, 0, l};

// Line connecting semi-circle's center point to left corner
L1 = newl;
Line(L1) = {C, P1};

// Segment connecting left corner to left candidates end-point
L2 = newl;
Circle(L2) = {P1, C, P2};

// Segment connecting candidates end-points
L3 = newl;
Circle(L3) = {P2, C, P3};

Transfinite Curve{L3} = n;

// Segment connecting right candidates end-point to right corner 
L4 = newl;
Circle(L4) = {P3, C, P4};

// Line connecting right corner to semi-circle's center point
L5 = newl;
Line(L5) = {P4, C};

Physical Curve("secondary_candidates") = {L3};
Physical Curve("secondary_bulk") = {L2, L4};
Physical Curve("secondary_fixed") = {L1, L5};

// Semi-circle
L = newl;
Curve Loop(L) = {L1, L2, L3, L4, L5};

S = newv;
Plane Surface(S) = {L};
Physical Surface("secondary") = {S};
