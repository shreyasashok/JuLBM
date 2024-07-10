function equilibrium(iPop::Int64, rho::Float64, u::AbstractArray{Float64}, uSqr::Float64, lat::Lattice)
    c_u::Float64 = 0.;
    for iD in 1:lat.d
        c_u += lat.c[iPop, iD] * u[iD];
    end
    return rho * lat.t[iPop] * (1. + lat.invCs2*c_u + lat.invCs2^2 * 0.5 * c_u^2 - lat.invCs2 * 0.5 * uSqr) - lat.t[iPop];
end
