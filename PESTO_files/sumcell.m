function out = sumcell(cellobj)
% sumcell sums the different elements of an 1-dimensional cell
% Cell should contain objects that can be added to each other (Points,
% functions, etc).
    out = cellobj{1};
    for i = 2:length(cellobj)
        out = out + cellobj{i};
    end
end

