function equilibriumBC!(data::AbstractArray{Float64}, lat::Lattice, omega::Float64, normOut::Vector{Int}) 
    rho = data[lat.rhoIndex];
    u = data[lat.uIndex:lat.uIndex+lat.d-1];
    uSqr = dot(u,u);
    data[1:lat.q] .= broadcast(x->equilibrium(x, rho, u, uSqr, D2Q9Lattice), 1:lat.q);
end