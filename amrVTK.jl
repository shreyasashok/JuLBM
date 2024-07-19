function amrVTK(grids::Vector{LBMData}, outFileName::String, it::Int, showFringe::Bool=true)
    baseSpacing::Float64 = 1.0;
    baseOrigin::Vector{Float64} = vec([0. 0.]);
    gridDict = Dict{Int64, Vector{Vector{Int64}}}();
    for grid in grids
        # println(grid.nx);
        if (grid.level == 0)
            baseSpacing = grid.dx;
            baseOrigin = vec([grid.plotX[1] grid.plotY[1]]);
        end
        levelVec = get!(gridDict, grid.level, []);

        plotXRange = 1:size(grid.plotX,1);
        plotYRange = 1:size(grid.plotY,1);
        xRange = 1:size(grid.plotX,1)-1;
        yRange = 1:size(grid.plotY,1)-1;
        fringeOffset = 0;
        if (!showFringe && grid.level > 0)
            plotXRange = 3:size(grid.plotX,1)-2;
            plotYRange = 3:size(grid.plotY,1)-2;
            xRange = 3:size(grid.plotX,1)-3;
            yRange = 3:size(grid.plotY,1)-3;
            fringeOffset = 2;
        end
        push!(levelVec, round.(Int, vec([(grid.plotX[1]-baseOrigin[1])/grid.dx+fringeOffset (grid.plotX[end]-baseOrigin[1])/grid.dx-1-fringeOffset (grid.plotY[1]-baseOrigin[2])/grid.dx+fringeOffset (grid.plotY[end]-baseOrigin[2])/grid.dx-1-fringeOffset 0 -1])));
        
        
        vtk_grid("vtk/$(outFileName)_level_$(grid.level)_grid$(length(levelVec)-1)_it$it", grid.plotX[plotXRange], grid.plotY[plotYRange]) do vtk
            vtk["rho"] = grid.data[D2Q9Lattice.rhoIndex, xRange, yRange];
            vtk["u"] = grid.data[D2Q9Lattice.uIndex+0, xRange, yRange];
            vtk["v"] = grid.data[D2Q9Lattice.uIndex+1, xRange, yRange];
            vtk["pop0"] = grid.data[1, xRange, yRange];
            vtk["pop1"] = grid.data[2, xRange, yRange];
            vtk["pop2"] = grid.data[3, xRange, yRange];
            vtk["pop3"] = grid.data[4, xRange, yRange];
            vtk["pop4"] = grid.data[5, xRange, yRange];
            vtk["pop5"] = grid.data[6, xRange, yRange];
            vtk["pop6"] = grid.data[7, xRange, yRange];
            vtk["pop7"] = grid.data[8, xRange, yRange];
            vtk["pop8"] = grid.data[9, xRange, yRange];
            vtk["pop2 plus pop6"] = grid.data[2, xRange, yRange] .+ grid.data[6, xRange, yRange];
        end
    end

    amrVTK = open("vtk/$(outFileName)_it$it.vthb", "w");
    println(amrVTK, "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">");
    println(amrVTK, "  <vtkOverlappingAMR origin=\"$(baseOrigin[1]) $(baseOrigin[2]) 0\" grid_description=\"XY\">");
    for level in gridDict
        levelSpacing = baseSpacing * 0.5^level.first;
        println(amrVTK, "    <Block level=\"$(level.first)\" spacing=\"$levelSpacing $levelSpacing $levelSpacing\">");
        for (igrid, grid) in enumerate(level.second)
            println(amrVTK, "      <DataSet index=\"$(igrid-1)\" amr_box=\"$(grid[1]) $(grid[2]) $(grid[3]) $(grid[4]) $(grid[5]) $(grid[6])\" file=\"$(outFileName)_level_$(level.first)_grid$(igrid-1)_it$it.vti\"/>");
        end
        println(amrVTK, "    </Block>");
    end
    println(amrVTK, "  </vtkOverlappingAMR>");
    println(amrVTK, "</VTKFile>");
    close(amrVTK);
end