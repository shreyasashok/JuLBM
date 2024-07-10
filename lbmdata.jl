struct LBMData

    nx::Int
    ny::Int
    data::Array{Float64}
    wrk::Array{Float64}
    x::Float64
    y::Float64
    dx::Float64
    plotX::AbstractRange{Float64}
    plotY::AbstractRange{Float64}
    cellX::AbstractRange{Float64} # cell centered data
    cellY::AbstractRange{Float64}
    

    function LBMData(nx::Int, ny::Int, lat::Lattice, x::Float64=0., y::Float64=0., dx::Float64=1.) 
        data = zeros(lat.dataSize, nx, ny)
        wrk = zeros(lat.d, nx, ny)
        plotX = (0:nx).*dx .+ x;
        plotY = (0:ny).*dx .+ y;
        cellX = (0:nx-1).*dx .+ (x + 0.5*dx)
        cellY = (0:ny-1).*dx .+ (y + 0.5*dx)
        new(nx,ny,data,wrk,x+0.5*dx,y+0.5*dx,dx,plotX,plotY,cellX,cellY)
    end

end