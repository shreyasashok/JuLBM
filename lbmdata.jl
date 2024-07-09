mutable struct LBMData

    const nx::Int
    const ny::Int
    data::Vector{Matrix{Float64}}
    const x::Float64
    const y::Float64
    const dx::Float64
    const plotX::AbstractRange{Float64}
    const plotY::AbstractRange{Float64}
    const cellX::AbstractRange{Float64} # cell centered data
    const cellY::AbstractRange{Float64}

    function LBMData(nx::Int, ny::Int, lat::Lattice, x::Float64=0., y::Float64=0., dx::Float64=1.) 
        data = Vector{Matrix{Float64}}(undef, lat.dataSize)
        for i in eachindex(data)
            data[i] = zeros(nx, ny)
        end
        plotX = (0:nx).*dx .+ x;
        plotY = (0:ny).*dx .+ y;
        cellX = (0:nx-1).*dx .+ (x + 0.5*dx)
        cellY = (0:ny-1).*dx .+ (y + 0.5*dx)
        new(nx,ny,data,x+0.5*dx,y+0.5*dx,dx,plotX,plotY,cellX,cellY)
    end

end