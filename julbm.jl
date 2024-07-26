using WriteVTK;
using Plots;

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

cellData = LBMData(5000,4,D2Q9Lattice, -2500., 0., 1.0);
push!(grids, cellData);
desiredViscosity = 10.0^-6;
tau = desiredViscosity*3 + 0.5;
push!(omegas, 1/tau);
iniShear!(cellData, D2Q9Lattice);

# collision! = bgkCollision!;
collision! = mrtCollision!;

nsteps = 90;
nsave = 1;

s1_even = 1.1;
s1_odd = 1.1;
s2_even = 1.1;
s2_odd = 1.1;

# s1_even = omegas[1];
# s1_odd = omegas[1];
# s2_even = omegas[1];
# s2_odd = omegas[1];

for i in 0:nsteps
    println("Step $i");
    if (mod(i,2) == 0)
        collision!(cellData, D2Q9Lattice, omegas[1], omegas[1], s1_even, s2_even)
    else
        collision!(cellData, D2Q9Lattice, omegas[1], omegas[1], s1_odd, s2_odd)
    end
    if (mod(i, nsave) == 0)
        amrVTK(grids, "vtkOut", i, true);
        print("Printing solution.\n")
    end
    stream!(cellData, D2Q9Lattice);
end

Uy = 0.1*sqrt(1.0/D2Q9Lattice.invCs2);
shearPlot = plot(cellData.cellX, cellData.data[D2Q9Lattice.uIndex+1, :, 2]/Uy)
xlabel!("x");
ylabel!("uy/Uy");
ylims!(-3E-3, 3E-3);
xlims!(0, 350);

# collision!(cellData, D2Q9Lattice, omegas[1])
collision!(cellData, D2Q9Lattice, omegas[1], omegas[1], s1_even, s2_even)
stream!(cellData, D2Q9Lattice);

plot!(cellData.cellX, cellData.data[D2Q9Lattice.uIndex+1, :, 2]/Uy)


print("Case complete.\n")
