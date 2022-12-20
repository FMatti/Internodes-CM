l = 0.3; // Target mesh size

r1 = 0.5; // Radius of the primary circle
r2 = 0.5; // Radius of the secondary circle
d = 0.05; // Initial penetration overlap of circles

a = 0.4; // Maximum expected radius of contact area
n = 20; // Number of mesh refinements for candidates


// ---> PRIMARY

// Opening angle for candidate nodes
phi = Asin(a/r1);

// Semi-circle center point
C = newp;
Point(C) = {0, -r1+d/2, 0, l};

// Semi-circle left corner point
P1 = newp;
Point(P1) = {r1, -r1+d/2, 0, l};

// Semi-circle candidate left end-point
P2 = newp;
Point(P2) = {r1*Sin(phi), r1*Cos(phi) -r1+d/2, 0, l};

// Circle candidate right end-point
P3 = newp;
Point(P3) = {-r1*Sin(phi), r1*Cos(phi) -r1+d/2, 0, l};

// Semi-circle right corner point
P4 = newp;
Point(P4) = {-r1, -r1+d/2, 0, l};

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

Physical Curve("primary_candidates") = {L3};
Physical Curve("primary_bulk") = {L2, L4};
Physical Curve("primary_fixed") = {L1, L5};

// Semi-circle
L = newl;
Curve Loop(L) = {L1, L2, L3, L4, L5};

S = newv;
Plane Surface(S) = {L};
Physical Surface("primary") = {S};


// ---> SECONDARY

// Opening angle for candidate nodes
phi = Asin(a/r2);

// Semi-circle center point
C = newp;
Point(C) = {0, r2-d/2, 0, l};

// Semi-circle left corner point
P1 = newp;
Point(P1) = {-r2, r2-d/2, 0, l};

// Semi-circle candidate left end-point
P2 = newp;
Point(P2) = {-r2*Sin(phi), -r2*Cos(phi) + r2-d/2, 0, l};

// Circle candidate right end-point
P3 = newp;
Point(P3) = {r2*Sin(phi), -r2*Cos(phi) + r2-d/2, 0, l};

// Semi-circle right corner point
P4 = newp;
Point(P4) = {r2, r2-d/2, 0, l};

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