// Gmsh project created on Mon Apr 15 18:48:03 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 2.3, 0.25, 1.45};
//+
Box(2) = {0.25, 0, 0, 1.8, 0.25, 1.2};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
