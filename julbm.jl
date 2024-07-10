using Plots; pythonplot();
using WriteVTK;

print("Hello world!\n")

include("lattice.jl");
include("latticeDefinitions.jl");
include("lbmdata.jl");
include("lbmhelpers.jl")
include("lbmops.jl");

cellData = LBMData(100,100,D2Q9Lattice, -50., -50., 1.0);
# iniEquilibrium(1.0, [0.1, 0.2], D2Q9Lattice, cellData);
iniVortex!(cellData, D2Q9Lattice);
# bgkCollision(D2Q9Lattice, cellData, 1.9);
# stream(D2Q9Lattice, cellData);

vtk_grid("test_it0", cellData.plotX, cellData.plotY) do vtk
    vtk["rho"] = cellData.data[D2Q9Lattice.rhoIndex];
    vtk["u"] = cellData.data[D2Q9Lattice.uIndex+0];
    vtk["v"] = cellData.data[D2Q9Lattice.uIndex+1];
end

for i in 1:5000
    bgkCollision!(cellData, D2Q9Lattice, 1.9);
    stream!(cellData, D2Q9Lattice);
    println(i);
    if (mod(i, 50) == 0)
        vtk_grid("test_it$i", cellData.plotX, cellData.plotY) do vtk
            vtk["rho"] = cellData.data[D2Q9Lattice.rhoIndex];
            vtk["u"] = cellData.data[D2Q9Lattice.uIndex+0];
            vtk["v"] = cellData.data[D2Q9Lattice.uIndex+1];
        end
    end
end