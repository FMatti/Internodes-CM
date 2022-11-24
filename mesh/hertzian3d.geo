l = 0.3;  // Target mesh size 
r = 0.75; // Radius of parabola at top
h = 0.3;  // Maximum height of candidate nodes
d = 0.05; // Penetration depth of parabola

// Discretization parameters
ny = 10;
nt = 5;
Geometry.ExtrudeSplinePoints = nt;
dy = r / (ny - 1);


// ---> PRIMARY
/*
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
*/

// ---> SECONDARY

// Generate parabola
ny = 2*r / l; // Number of nodes along parabola
dx = r / (ny - 1); // Spacing between nodes in x-direction
j = 0;
For i In {0:ny-1}

    // Generate new point on parabola
    x = i*dx;
    y0 = 1 - d;
    y = y0 + x*x;

    p = newp;
    Point(p) = {x, y, 0, l};
    P[i] = p;

    // Connect points with lines
    If (i > 0)
        l = newl;
        Line(l) = {P[i-1], P[i]};
        L[i-1] = l;
        
        // Add to candidate list if below height h
        If (y < y0 + h)
            C[j] = l;
            j += 1;
        EndIf
    EndIf
EndFor

// Close off the top
p = newp;
Point(p) = {0, y, 0, l};
l = newl;
Line(l) = {P[ny-1], p};
Transfinite Curve{l} = ny;
L[ny-1] = l;

For i In {0:ny-1}
    extrusion[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{L[i]}; };
    extrusion2[] = Extrude{{0, 0, 0}, {0, -1, 0}, {0, 0, 0}, Pi}{ Line{extrusion[0]}; };
    S[2*i] = extrusion[1];
    S[2*i + 1] = extrusion2[1];
    If (i == 0)
        Physical Surface("secondary_candidate") = extrusion[1];
        Physical Surface("secondary_candidate") += extrusion2[1];
    ElseIf (i == ny-1)
        Physical Surface("secondary_fixed") = extrusion[1];
        Physical Surface("secondary_fixed") += extrusion2[1];
    ElseIf (i > 0)
        Physical Surface("secondary_candidate") += extrusion[1];
        Physical Surface("secondary_candidate") += extrusion2[1];
    EndIf
EndFor

Surface Loop(2) = {S};
Volume(2) = {2};
Physical Volume("secondary") = {2};