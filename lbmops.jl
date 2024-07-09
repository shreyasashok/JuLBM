using LinearAlgebra;

function iniEquilibrium(rho::Float64, u::Vector{Float64}, lat::Lattice, cellData::LBMData)
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

function bgkCollision(lat::Lattice, cellData::LBMData, omega::Float64)
    
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

function stream(lat::Lattice, cellData::LBMData)
    for iPop in 1:lat.q
        circshift(cellData.data[iPop], lat.c[iPop,:])
    end
end