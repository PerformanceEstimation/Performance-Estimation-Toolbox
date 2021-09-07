function [G,F]=LocalOracles(Fi,X)
% LocalOracles call an oracle for each objective function contained in the cell Fi
% to obtain the value and the gradient of the functions for a different point contained in cell X.
% The different functions and points typically corresponds to the
% differents local functions and variables of agents in decentralized optimization. 
% Input :   - Fi is a cell (1xN) of objective functions (from FunctionClass)
%           - X is a cell (1xN) of Points or Evaluable object. 
%           It should have the same size as Fi.
% Output :  G and F are two cells (1xN) containing the function value and
%           the gradient value of each objective function Fi, for a different point
%           from cell X.
%           Gi(Xi) and Fi(Xi)
    assert(length(Fi)==length(X),'Fi and X should be two cells of the same size.');
    N = length(X);
    G = cell(N,1); F = cell(N,1);
    for i=1:N
        [G{i}, F{i}] = Fi{i}.oracle(X{i});
    end
end