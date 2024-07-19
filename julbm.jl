using WriteVTK;

print("Starting JuLBM\n")

include("lattice.jl");
include("latticeDefinitions.jl");
include("lbmdata.jl");
include("lbmhelpers.jl")
include("lbmops.jl");
include("amrVTK.jl");
include("fringeInterp.jl")
include("bc.jl")

rm("vtk/", force=true, recursive=true);
mkpath("vtk/");

grids::Vector{LBMData} = [];
omegas::Vector{Float64} = [];

lbmData = LBMData(250,90,D2Q9Lattice, -50.5, -50., 1.0);
lbmDataFine = LBMData(70,80,D2Q9Lattice, 50.5, -30., 0.5, 1);
push!(grids, lbmData);
push!(grids, lbmDataFine);
push!(omegas, 1.9);
push!(omegas, omegaForLevel(1.9, 1));
# iniEquilibrium!(lbmData, D2Q9Lattice, 1.0, [0.1, 0.2]);
# iniEquilibrium!(lbmDataFine, D2Q9Lattice, 1.0, [-0.1, -0.2]);
iniVortex!(lbmDataFine, D2Q9Lattice);
iniVortex!(lbmData, D2Q9Lattice);

lbmData.data[1, 125, 45] = NaN;

# collision! = bgkCollision!;
collision! = mrtCollision!;

bcs = [equilibriumBC!, #xMin
       equilibriumBC!, #xMax
       equilibriumBC!, #yMin
       equilibriumBC!, #yMax
       equilibriumBC!, #xMin yMin
       equilibriumBC!, #xMax yMin
       equilibriumBC!, #xMin yMax
       equilibriumBC!]; #xMax yMax

nsteps = 1500;
nsave = nsteps/30;

for i in 0:nsteps
    println("Step $i");
    collision!(lbmData, D2Q9Lattice, omegas[1])
    collision!(lbmDataFine, D2Q9Lattice, omegas[2]);
    fineToCoarse!(lbmData, lbmDataFine);
    coarseToFine!(lbmDataFine, lbmData);
    if (mod(i, nsave) == 0)
        amrVTK(grids, "vtkOut", i, true);
        print("Printing solution.\n")
    end
    stream!(lbmData, D2Q9Lattice);
    # applyBCs!(lbmData, D2Q9Lattice, omegas[1], bcs);
    stream!(lbmDataFine, D2Q9Lattice);
    if (mod(i, nsave) == 0)
        amrVTK(grids, "vtkOut", i+1, true);
        print("Printing solution.\n")
    end
    collision!(lbmDataFine, D2Q9Lattice, omegas[2]);
    if (mod(i, nsave) == 0)
        amrVTK(grids, "vtkOut", i+2, true);
        print("Printing solution.\n")
    end
    stream!(lbmDataFine, D2Q9Lattice);
    if (mod(i, nsave) == 0)
        amrVTK(grids, "vtkOut", i+3, true);
        print("Printing solution.\n")
    end
end

print("Case complete.\n")
