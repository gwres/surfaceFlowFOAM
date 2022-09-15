lc =  1.0;

//Points
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {1.0, 400.0, 0.0, lc};
Point(4) = {0.0, 400.0, 0.0, lc};

//Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Line {1, 3} = 6 Using Progression 1.0;
Transfinite Line {2, 4} = 101 Using Progression 1.0;

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 1.0} {
	Surface{1};
	Layers{1};
	Recombine;
	}
//+
Physical Surface("inlet") = {13};
//+
Physical Surface("noFlow") = {25, 17};
//+
Physical Surface("outlet") = {21};
//+
Physical Volume("domain") = {1};
