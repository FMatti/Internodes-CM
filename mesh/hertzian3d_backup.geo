l = 0.2;  // Mesh size
r = 1.0;  // Radius of paraboloid at highest point
y0 = 0.9;  // Position of tip of paraboloid

// Discretization parameters
ny = 10;
nt = 10;
dy = r / (ny - 1);

// ---> PRIMARY

// Nodes
Point(1) = {-0.5, 0, -0.5, l};
Point(2) = {0.5, 0, -0.5, l};
Point(3) = {0.5, 1, -0.5, l};
Point(4) = {-0.5, 1, -0.5, l};
Point(5) = {-0.5, 0, 0.5, l};
Point(6) = {0.5, 0, 0.5, l};
Point(7) = {0.5, 1, 0.5, l};
Point(8) = {-0.5, 1, 0.5, l};

// Lines
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

// Surfaces
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
Physical Surface("primary_candidate") = {6};

// Volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Volume("primary") = {1};


// ---> SECONDARY

// Origin point
P0 = newp;
P[0] = P0;
Point(P0) = {0, y0, 0, l};

For p In {1:ny}
    x = p*dy;
    If (p < ny)
        // Parabolic function
        y = y0 + x*x;
    Else
        // Set final point to 0 to close off surface
        x = 0;
    EndIf

    // Generate new point on surface
    Pp = newp;
    Point(Pp) = {x, y, 0, l};
    P[p] = Pp;

    // Connect points with line
    L = newl;
    Line(L) = {P[p-1], P[p]};

    // Revolution-extrude line segment twice (two half-turns)
    extrusion[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{L}; };
    extrusion2[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{extrusion[0]}; };

    If (p == 1)
        Physical Surface("secondary_candidate") = extrusion[1];
        Physical Surface("secondary_candidate") += extrusion2[1];
    ElseIf (p == ny)
        Physical Surface("secondary_fixed") = extrusion[1];
        Physical Surface("secondary_fixed") += extrusion2[1];
    ElseIf (p > 1)
        Physical Surface("secondary_candidate") += extrusion[1];
        Physical Surface("secondary_candidate") += extrusion2[1];
    EndIf
EndFor

Surface Loop(2) = {21, 30, 39, 48, 57, 66, 75, 84, 93, 17, 26, 35, 44, 53, 62, 71, 80, 89, 98, 102};
Volume(2) = {2};
Physical Volume("secondary") = {2};
