function out = sumcell(cellobj)
% sumcell sums the different elements of an cell
% Cell should contain objects that can be added to each other (Points,
% functions, etc).
    cell_size = size(cellobj);
    out = cellobj{1};
    for i=2:cell_size(1)*cell_size(2)
            out = out + cellobj{i};
    end
end