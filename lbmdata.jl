mutable struct LBMData

    const nx::Int
    const ny::Int
    data::Vector{Matrix{Float64}}

    function LBMData(nx::Int, ny::Int, lat::Lattice) 
        data = Vector{Matrix{Float64}}(undef, lat.dataSize)
        for i in eachindex(data)
            data[i] = zeros(nx, ny)
        end
        new(nx,ny,data)
    end
end