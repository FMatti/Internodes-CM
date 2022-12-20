l = 0.3; // Target mesh size 

r1 = 0.5; // Radius of primary sphere (bottom)
r2 = 0.5; // Radius of secondary sphere (top)
d = 0.05; // Initial penetration overlap of sphere

a = 0.4; // Maximum expected radius of contact area
n = 10; // Number of mesh refinements for candidates


// ---> PRIMARY

Geometry.ExtrudeSplinePoints = 10;

// Opening angle for candidate nodes
phi = Asin(a/r1); 

// Sphere center point
C = newp;
Point(C) = {0, -r1+d/2, 0, l};

// Sphere top point
P0 = newp;
Point(P0) = {0, d/2, 0, l};

// Sphere candidates end-point
P1 = newp;
Point(P1) = {-r1*Sin(phi), r1*Cos(phi) -r1+d/2, 0, l};

// Sphere end-point
P2 = newp;
Point(P2) = {-r1, -r1+d/2, 0, l};

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
E1a[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, -Pi}{ Line{L1}; };
E1b[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, -Pi}{ Line{E1a[0]}; };

// Extrude bulk surface
E2a[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, -Pi}{ Line{L2}; };
E2b[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, -Pi}{ Line{E2a[0]}; };

// Extrude top surface
E3a[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, -Pi}{ Line{L3}; };
E3b[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, -Pi}{ Line{E3a[0]}; };

Transfinite Curve(L1) = n;
Transfinite Curve(E1a[0]) = n;
Transfinite Curve(E1a[2]) = 2*n;
Transfinite Curve(E1b[2]) = 2*n;

Physical Curve("primary_transfinite") = {L1, E1a[0], E1a[2], E1b[2]};

Physical Surface("primary_candidates") = {E1a[1], E1b[1]};
Physical Surface("primary_bulk") = {E2a[1], E2b[1]};
Physical Surface("primary_fixed") = {E3a[1], E3b[1]};

// Semi-sphere
S = news;
Surface Loop(S) = {E1a[1], E1b[1], E2a[1], E2b[1], E3a[1], E3b[1]};

V = newv;
Volume(V) = {S};
Physical Volume("primary") = {V};

// ---> SECONDARY

Geometry.ExtrudeSplinePoints = 10;

// Opening angle for candidate nodes
phi = Asin(a/r2); 

// Sphere center point
C = newp;
Point(C) = {0, r2-d/2, 0, l};

// Sphere bottom point
P0 = newp;
Point(P0) = {0, -d/2, 0, l};

// Sphere candidates end-point
P1 = newp;
Point(P1) = {r2*Sin(phi), -r2*Cos(phi) + r2-d/2, 0, l};

// Sphere end-point
P2 = newp;
Point(P2) = {r2, r2-d/2, 0, l};

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

Transfinite Curve(L1) = n;
Transfinite Curve(E1a[0]) = n;
Transfinite Curve(E1a[2]) = 2*n;
Transfinite Curve(E1b[2]) = 2*n;

Physical Curve("secondary_transfinite") = {L1, E1a[0], E1a[2], E1b[2]};

Physical Surface("secondary_candidates") = {E1a[1], E1b[1]};
Physical Surface("secondary_bulk") = {E2a[1], E2b[1]};
Physical Surface("secondary_fixed") = {E3a[1], E3b[1]};

// Semi-sphere
S = news;
Surface Loop(S) = {E1a[1], E1b[1], E2a[1], E2b[1], E3a[1], E3b[1]};

V = newv;
Volume(V) = {S};
Physical Volume("secondary") = {V};
