struct Lattice
    d::Int
    q::Int
    rhoIndex::Int
    uIndex::Int 
    dataSize::Int
    c::Matrix{Int64}
    t::Vector{Float64}
    invCs2::Float64
    M::AbstractMatrix{Float64}
    Minv::AbstractMatrix{Float64}
end
