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

cellData = LBMData(250,90,D2Q9Lattice, -50.5, -50., 1.0);
cellDataFine = LBMData(70,80,D2Q9Lattice, 50.5, -30., 0.5, 1);
push!(grids, cellData);
push!(grids, cellDataFine);
push!(omegas, 1.9);
push!(omegas, omegaForLevel(1.9, 1));
# iniEquilibrium!(cellData, D2Q9Lattice, 1.0, [0.1, 0.2]);
# iniEquilibrium!(cellDataFine, D2Q9Lattice, 1.0, [-0.1, -0.2]);
iniVortex!(cellDataFine, D2Q9Lattice);
iniVortex!(cellData, D2Q9Lattice);

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

nsteps = 5000;
nsave = nsteps/10;

for i in 0:nsteps
    println("Step $i");
    collision!(cellData, D2Q9Lattice, omegas[1])
    collision!(cellDataFine, D2Q9Lattice, omegas[2]);
    coarseToFine!(cellDataFine, cellData);
    fineToCoarse!(cellData, cellDataFine);
    if (mod(i, nsave) == 0)
        amrVTK(grids, "vtkOut", i, true);
        print("Printing solution.\n")
    end
    stream!(cellData, D2Q9Lattice);
    applyBCs!(cellData, D2Q9Lattice, omegas[1], bcs);
    stream!(cellDataFine, D2Q9Lattice);
    collision!(cellDataFine, D2Q9Lattice, omegas[2]);
    stream!(cellDataFine, D2Q9Lattice);
end

print("Case complete.\n")
