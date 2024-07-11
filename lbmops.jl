using LinearAlgebra;

function iniEquilibrium!(cellData::LBMData, lat::Lattice, rho::Float64, u::Vector{Float64})
    uSqr::Float64 = dot(u,u);
    for idx in CartesianIndices(cellData.data[1,:,:])
        for iPop in 1:lat.q
            cellData.data[iPop,idx] = equilibrium(iPop, rho, u, uSqr, lat);
        end
        cellData.data[lat.rhoIndex, idx] = rho;
        for iD in 1:lat.d
            cellData.data[lat.uIndex-1+iD, idx] = u[iD];
        end
    end
end

function iniVortex!(cellData::LBMData, lat::Lattice)
    for idx in CartesianIndices(cellData.data[1,:,:])
        i, j = Tuple.(idx);
        coords = [cellData.cellX[i], cellData.cellY[j]];
        radSqr = dot(coords, coords);
        vortexRadius::Float64 = 5.0;
        ee = exp(0.5*(-radSqr/vortexRadius^2));
        strength::Float64 = 0.01;
        u = 0.05 - strength*ee*coords[2];
        v = strength*ee*coords[1];

        rho = 1.0;
        uVec = [u,v];
        uSqr = dot(uVec,uVec);

        for iPop in 1:lat.q
            cellData.data[iPop,idx] = equilibrium(iPop, rho, uVec, uSqr, lat);
        end
        cellData.data[lat.rhoIndex,idx] = rho;
        for iD in 1:lat.d
            cellData.data[lat.uIndex-1+iD, idx] = uVec[iD];
        end

    end
end

function bgkCollision!(cellData::LBMData, lat::Lattice, omega::Float64) 
    vectorizedBGK(cellData.data, cellData.wrk, lat, omega);
end

function vectorizedBGK(data, wrk, lat::Lattice, omega::Float64) 
    for idx in CartesianIndices(@view(data[1,:,:]))
        bgkCollisionOperation!(@view(data[:,idx]), @view(wrk[:,idx]), lat, omega);
    end
end


function bgkCollisionOperation!(data::AbstractArray{Float64}, uWrk::AbstractArray{Float64}, lat::Lattice, omega::Float64) 
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

function stream!(cellData::LBMData, lat::Lattice)
    for iPop in 1:lat.q
        cellData.data[iPop,:,:] = circshift(cellData.data[iPop,:,:], lat.c[iPop,:])
    end
end