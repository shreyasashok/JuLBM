using WriteVTK;

print("Hello world!\n")

include("lattice.jl");
include("latticeDefinitions.jl");
include("lbmdata.jl");
include("lbmhelpers.jl")
include("lbmops.jl");
include("amrVTK.jl");

rm("vtk/", force=true, recursive=true);
mkpath("vtk/");

grids::Vector{LBMData} = [];

cellData = LBMData(100,90,D2Q9Lattice, -50.5, -50., 1.0);
cellDataFine = LBMData(70,80,D2Q9Lattice, -40.5, -30., 0.5, 1);
push!(grids, cellData);
push!(grids, cellDataFine);
# iniEquilibrium(1.0, [0.1, 0.2], D2Q9Lattice, cellData);
iniVortex!(cellData, D2Q9Lattice);
iniVortex!(cellDataFine, D2Q9Lattice);

for i in 1:5000
    bgkCollision!(cellData, D2Q9Lattice, 1.9);
    stream!(cellData, D2Q9Lattice);
    println(i);
    if (mod(i, 50) == 0)
        amrVTK(grids, "vtkOut", i, false);
    end
end