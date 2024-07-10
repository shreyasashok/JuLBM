function testAlloc(u::Vector{Float64})
    fill!(u, 0.0);
    return nothing;
end

test = vec(ones(10000,1));

include("lattice.jl");
include("latticeDefinitions.jl");
include("lbmdata.jl");
include("lbmhelpers.jl")
include("lbmops.jl");

function bgkCollisionOperationTest!(data::AbstractArray{Float64}, uWrk::AbstractArray{Float64}, lat::Lattice, omega::Float64) 
    rho::Float64 = 1.0;
    fill!(uWrk, 0.0);

    for iPop in 1:lat.q
        rho += data[iPop];
        for iD in 1:lat.d
            uWrk[iD] += data[iPop] * lat.c[iPop, iD];
        end
    end
    uWrk ./= rho;
    uSqr = dot(uWrk,uWrk);

    data[lat.rhoIndex] = rho;
    for iD in 1:lat.d
        data[lat.uIndex-1+iD] = uWrk[iD];
    end
    
    for iPop in 1:lat.q
        data[iPop] = (data[iPop] * (1.0 - omega)) + omega*equilibrium(iPop, rho, uWrk, uSqr, lat);
    end
    return nothing;
end

function vectorizedBGK(data, wrk, lat::Lattice, omega::Float64) 
    for idx in CartesianIndices(@view(data[1,:,:]))
        bgkCollisionOperationTest!(@view(data[:,idx]), @view(wrk[:,idx]), lat, omega);
    end
end

testPops = vec(rand(12, 1));
testWrk = vec(rand(2,1));

cellData = LBMData(100,100,D2Q9Lattice, -50., -50., 1.0);
wrk = zeros(2,cellData.nx, cellData.ny);

