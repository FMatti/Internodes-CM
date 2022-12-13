l = 0.1;  // Target mesh size 
r = 0.5; // Radius of parabola at top
h = 0.2;  // Maximum height of candidate nodes
d = 0.05; // Penetration depth of parabola


// ---> PRIMARY

// Nodes
Point(1) = {-0.5, 0, 0, l};
Point(2) = {0.5, 0, 0, l};
Point(3) = {0.5, 0.5, 0, l};
Point(4) = {-0.5, 0.5, 0, l};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Physical Curve("primary_fixed") = {1};
Physical Curve("primary_candidate") = {3};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("primary") = {1};


// ---> SECONDARY

// Generate parabola
ny = 2*r / l; // Number of nodes along parabola
dx = (2*r) / (ny - 1); // Spacing between nodes in x-direction
j = 0;
For i In {0:ny-1}

    // Generate new point on parabola
    x = i*dx - r;
    y0 = 0.5 - d;
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
l = newl;
Line(l) = {P[ny-1], P[0]};
Transfinite Curve{l} = ny;
L[ny-1] = l;

Physical Curve("secondary_candidate") = {C[1]:C[j-1]};
Physical Curve("secondary_fixed") = {l};

// Create surface
Curve Loop(2) = {L[0]:L[ny-1]};
Plane Surface(2) = {2};
Physical Surface("secondary") = {2};
