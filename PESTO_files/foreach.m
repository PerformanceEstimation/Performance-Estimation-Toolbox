function out=foreach(func, varargin)
% foreach apply the same function 'func' for each elements of the inputs cells contained in varargin
    out = cellfun(func, varargin{:},'UniformOutput',false);
end