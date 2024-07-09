using Plots; pythonplot();

print("Hello world!\n")

include("lattice.jl");
include("latticeDefinitions.jl");
include("lbmdata.jl");
include("lbmhelpers.jl")
include("lbmops.jl");

cellData = LBMData(5,6,D2Q9Lattice);
iniEquilibrium(1.0, [0.1, 0.2], D2Q9Lattice, cellData);
bgkCollision(D2Q9Lattice, cellData, 1.9);
stream(D2Q9Lattice, cellData);