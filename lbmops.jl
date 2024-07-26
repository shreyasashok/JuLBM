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

function iniShear!(cellData::LBMData, lat::Lattice)
    cs = sqrt(1.0/lat.invCs2);
    Ux = 0.1*cs;
    Uy = 0.1*cs;
    rho = 1.0;
    xc = 150;
    Rc = 5;
    for idx in CartesianIndices(cellData.data[1,:,:])
        x = cellData.cellX[idx[1]];
        ux = Ux;
        uy = Uy * exp(-((x-xc)^2)/(2*Rc^2));
        u = [ux,uy];
        uSqr = dot(u,u);       

        for iPop in 1:lat.q
            cellData.data[iPop,idx] = equilibrium(iPop, rho, u, uSqr, lat);
        end
        cellData.data[lat.rhoIndex, idx] = rho;
        for iD in 1:lat.d
            cellData.data[lat.uIndex-1+iD, idx] = u[iD];
        end
    end
end

function bgkCollision!(cellData::LBMData, lat::Lattice, omega::Float64) 
    vectorizedBGK!(cellData.data, cellData.wrk, lat, omega);
end

function mrtCollision!(cellData::LBMData, lat::Lattice, omega::Float64, bulk::Float64=omega, s1::Float64=1.1, s2::Float64=1.1) 
    vectorizedMRT!(cellData.data, cellData.wrk, cellData.mWrk, lat, omega, lat.M, lat.Minv, bulk, s1, s2);
end

function vectorizedBGK!(data, wrk, lat::Lattice, omega::Float64) 
    for idx in CartesianIndices(@view(data[1,:,:]))
        bgkCollisionOperation!(@view(data[:,idx]), @view(wrk[:,idx]), lat, omega);
    end
end

function vectorizedMRT!(data, uWrk, mWrk, lat::Lattice, omega::Float64, M::AbstractArray{Float64}, Minv::AbstractArray{Float64}, bulk::Float64, s1::Float64, s2::Float64) 
    for idx in CartesianIndices(@view(data[1,:,:]))
        # if mod(mod(idx[1], 2) + mod(idx[2],2), 2) == 0
            mrtCollisionOperation!(@view(data[:,idx]), @view(uWrk[:,idx]), @view(mWrk[:,idx]), lat, M, Minv, omega, bulk, s1, s2);
        # else
            # mrtCollisionOperation!(@view(data[:,idx]), @view(uWrk[:,idx]), @view(mWrk[:,idx]), lat, M, Minv, omega, bulk, 2.0-s1, 2.0-s2);
        # end
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

function mrtCollisionOperation!(data::AbstractArray{Float64}, uWrk::AbstractArray{Float64}, mWrk::AbstractArray{Float64}, lat::Lattice, M::AbstractArray{Float64}, Minv::AbstractArray{Float64}, omega::Float64, bulk::Float64, s1::Float64, s2::Float64) 
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

    mWrk .= @view(data[1:lat.q]) .- broadcast(x->equilibrium(x, rho, uWrk, uSqr, D2Q9Lattice), 1:lat.q); #fNeq
    mWrk .= M*mWrk; #mNeq
    mWrk = Array(mWrk); #convert from SVector to normal vector
    mWrk[4:4] .*= -bulk; #relax (bulk viscosity)
    mWrk[5:5] .*= -s1; #relax (free parameter)
    mWrk[6:7] .*= -s2; #relax (free parameter)
    mWrk[8:9] .*= -omega; #relax (shear viscosity)

    data[lat.rhoIndex] = rho;
    for iD in 1:lat.d
        data[lat.uIndex-1+iD] = uWrk[iD];
    end

    data[1:lat.q] .+= Minv*mWrk;
    return nothing;
end

function stream!(cellData::LBMData, lat::Lattice)
    vectorizedStream!(cellData.data, lat);
end

function vectorizedStream!(data::AbstractArray{Float64}, lat::Lattice)
    for iPop in 1:lat.q
        data[iPop,:,:] = circshift(@view(data[iPop,:,:]), lat.c[iPop,:])
    end
end

#order of bcs: xMin, xMax, yMin, yMax, xMinyMin, xMaxyMin, xMinYMax, xMaxYMax
function applyBCs!(cellData::LBMData, lat::Lattice, omega::Float64, bcfs!)
    applyBC!(cellData, lat, omega, bcfs![1], [-1,  0]);
    applyBC!(cellData, lat, omega, bcfs![2], [ 1,  0]);
    applyBC!(cellData, lat, omega, bcfs![3], [ 0, -1]);
    applyBC!(cellData, lat, omega, bcfs![4], [ 0,  1]);
    applyBC!(cellData, lat, omega, bcfs![5], [-1, -1]);
    applyBC!(cellData, lat, omega, bcfs![6], [ 1, -1]);
    applyBC!(cellData, lat, omega, bcfs![7], [-1,  1]);
    applyBC!(cellData, lat, omega, bcfs![8], [ 1,  1]);
end

function applyBC!(cellData::LBMData, lat::Lattice, omega::Float64, bcf!::Function, normOut::Vector{Int})
    xRange = normOut[1] == 0 ? (2:cellData.nx-1) : normOut[1] == (-1:-1) ? (1:1) : (cellData.nx:cellData.nx);
    yRange = normOut[2] == 0 ? (2:cellData.ny-1) : normOut[2] == (-1:-1) ? (1:1) : (cellData.ny:cellData.ny);

    vectorizedBC!(cellData.data, xRange, yRange, lat, omega, bcf!, normOut);
end

function vectorizedBC!(data, xRange, yRange, lat::Lattice, omega::Float64, bcf!::Function, normOut::Vector{Int})
    for idx in CartesianIndices(@view(data[1,xRange,yRange]))
        bcf!(@view(data[:,idx]), lat, omega, normOut);
    end
end