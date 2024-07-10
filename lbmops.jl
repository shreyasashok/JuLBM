using LinearAlgebra;

function iniEquilibrium!(cellData::LBMData, lat::Lattice, rho::Float64, u::Vector{Float64})
    uSqr::Float64 = dot(u,u);
    for idx in CartesianIndices(cellData.data[1])
        for iPop in 1:lat.q
            cellData.data[iPop][idx] = equilibrium(iPop, rho, u, uSqr, lat);
        end
        cellData.data[lat.rhoIndex][idx] = rho;
        for iD in 1:lat.d
            cellData.data[lat.uIndex-1+iD][idx] = u[iD];
        end
    end
end

function iniVortex!(cellData::LBMData, lat::Lattice)
    for idx in CartesianIndices(cellData.data[1])
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
            cellData.data[iPop][idx] = equilibrium(iPop, rho, uVec, uSqr, lat);
        end
        cellData.data[lat.rhoIndex][idx] = rho;
        for iD in 1:lat.d
            cellData.data[lat.uIndex-1+iD][idx] = uVec[iD];
        end

    end
end

function bgkCollision!(cellData::LBMData, lat::Lattice, omega::Float64)
    
    for idx in CartesianIndices(cellData.data[1])
        rho::Float64 = 1.0;
        u::Vector{Float64} = zeros(lat.q);
        for iPop in 1:lat.q
            rho += cellData.data[iPop][idx];
            for iD in 1:lat.d
                u[iD] += cellData.data[iPop][idx] * lat.c[iPop, iD];
            end
        end
        u /= rho;
        uSqr = dot(u,u);

        cellData.data[lat.rhoIndex][idx] = rho;
        for iD in 1:lat.d
            cellData.data[lat.uIndex-1+iD][idx] = u[iD];
        end
        
        for iPop in 1:lat.q
            cellData.data[iPop][idx] = (cellData.data[iPop][idx] * (1.0 - omega)) + omega*equilibrium(iPop, rho, u, uSqr, lat);
        end        
    end
end

function stream!(cellData::LBMData, lat::Lattice)
    for iPop in 1:lat.q
        cellData.data[iPop] = circshift(cellData.data[iPop], lat.c[iPop,:])
    end
end