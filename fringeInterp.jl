function getInterpolatedValuesOperator(data::AbstractArray{Float64}, points::AbstractArray{Float64}, interpGrid::RectangleGrid)
    return permutedims(cat(map(x -> map(y -> interpolate(interpGrid, x, y), eachslice(points, dims=(2,3))), eachslice(data, dims=1))..., dims=3), (3,1,2));
end

function getInterpolatedValues(cellData::LBMData, points::AbstractArray{Float64})
    return getInterpolatedValuesOperator(cellData.data, points, cellData.interpGrid);
end

function fineToCoarse!(cellDataCoarse::LBMData, cellDataFine::LBMData) 

    startXCoord = 0.5 * (cellDataFine.cellX[3] + cellDataFine.cellX[4]);
    startYCoord = 0.5 * (cellDataFine.cellY[3] + cellDataFine.cellY[4]);
    endXCoord = 0.5 * (cellDataFine.cellX[end-3] + cellDataFine.cellX[end-2]);
    endYCoord = 0.5 * (cellDataFine.cellY[end-3] + cellDataFine.cellY[end-2]);

    startXIndex = findfirst(abs.(cellData.cellX .- startXCoord) .< 1E-7)
    startYIndex = findfirst(abs.(cellData.cellY .- startYCoord) .< 1E-7)
    endXIndex = findfirst(abs.(cellData.cellX .- endXCoord) .< 1E-7)
    endYIndex = findfirst(abs.(cellData.cellY .- endYCoord) .< 1E-7)

    #Bottom edge
    cellDataCoarse.data[:, startXIndex:endXIndex, startYIndex] = getInterpolatedValues(cellDataFine, permutedims(cat([x for x=cellDataCoarse.cellX[startXIndex:endXIndex], y=cellDataCoarse.cellY[startYIndex]], 
                                                                                                                     [y for x=cellDataCoarse.cellX[startXIndex:endXIndex], y=cellDataCoarse.cellY[startYIndex]], dims=3),(3,1,2)));

    #Top edge
    cellDataCoarse.data[:, startXIndex:endXIndex, endYIndex] = getInterpolatedValues(cellDataFine, permutedims(cat([x for x=cellDataCoarse.cellX[startXIndex:endXIndex], y=cellDataCoarse.cellY[endYIndex]], 
                                                                                                                   [y for x=cellDataCoarse.cellX[startXIndex:endXIndex], y=cellDataCoarse.cellY[endYIndex]], dims=3),(3,1,2)));

    #Left edge minus top and bottom
    cellDataCoarse.data[:, startXIndex, startYIndex+1:endYIndex-1] = getInterpolatedValues(cellDataFine, permutedims(cat([x for x=cellDataCoarse.cellX[startXIndex], y=cellDataCoarse.cellY[startYIndex+1:endYIndex-1]], 
                                                                                                                         [y for x=cellDataCoarse.cellX[startXIndex], y=cellDataCoarse.cellY[startYIndex+1:endYIndex-1]], dims=3),(3,1,2)));

    #Right edge minus top and bottom
    cellDataCoarse.data[:, endXIndex, startYIndex+1:endYIndex-1] = getInterpolatedValues(cellDataFine, permutedims(cat([x for x=cellDataCoarse.cellX[endXIndex], y=cellDataCoarse.cellY[startYIndex+1:endYIndex-1]], 
                                                                                                                       [y for x=cellDataCoarse.cellX[endXIndex], y=cellDataCoarse.cellY[startYIndex+1:endYIndex-1]], dims=3),(3,1,2)));
end

function coarseToFine!(cellDataFine::LBMData, cellDataCoarse::LBMData) 

    #Bottom edge
    cellDataFine.data[:, :, 1:2] = getInterpolatedValues(cellDataCoarse, permutedims(cat([x for x=cellDataFine.cellX, y=cellDataFine.cellY[1:2]], 
                                                                                         [y for x=cellDataFine.cellX, y=cellDataFine.cellY[1:2]], dims=3),(3,1,2)));

    #Top edge
    cellDataFine.data[:, :, end-1:end] = getInterpolatedValues(cellDataCoarse, permutedims(cat([x for x=cellDataFine.cellX, y=cellDataFine.cellY[end-1:end]], 
                                                                                               [y for x=cellDataFine.cellX, y=cellDataFine.cellY[end-1:end]], dims=3),(3,1,2)));

    #Left edge minus top and bottom
    cellDataFine.data[:, 1:2, 3:end-2] = getInterpolatedValues(cellDataCoarse, permutedims(cat([x for x=cellDataFine.cellX[1:2], y=cellDataFine.cellY[3:end-2]], 
                                                                                               [y for x=cellDataFine.cellX[1:2], y=cellDataFine.cellY[3:end-2]], dims=3),(3,1,2)));

    #Right edge minus top and bottom
    cellDataFine.data[:, end-1:end, 3:end-2] = getInterpolatedValues(cellDataCoarse, permutedims(cat([x for x=cellDataFine.cellX[end-1:end], y=cellDataFine.cellY[3:end-2]], 
                                                                                                     [y for x=cellDataFine.cellX[end-1:end], y=cellDataFine.cellY[3:end-2]], dims=3),(3,1,2)));
end