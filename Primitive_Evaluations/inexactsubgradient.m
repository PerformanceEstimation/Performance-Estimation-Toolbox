function [ d ] = inexactsubgradient( x,f,eps,abs )
% [ d ] = inexactsubgradient( x,f,eps,abs )
%
% This routine allows to evaluate an inexact (sub)gradient. 
%
% Inputs:
%   - x is a point in the decision space (not a function value)
%   - f is a function
%   - eps is the required accuracy (either relative or absolute, see above)
%   - abs defines the mode (absolute or relative inaccuracy, see above)
%
% Relative or absolute accuracy; the evaluation can be performed in two
% modes:
%   - relative inaccuracy: (default)
%       ||d-g||<=eps*||g|| where g is a (sub)gradient of f at x
%   - absolute inaccuracy: (set abs=1) 
%       ||d-g||<=eps where g is a (sub)gradient of f at x
%
% Output: inexact (sub)gradient d
%
%
% Example 1; evaluate an inexact subgradient (relative accuracy) of a
% Function F at some point x0
%
% >> P=pet(); param.L=1;param.mu=0.2; eps=.1;
% >> F=P.AddObjective('SmoothStronglyConvex',param);
% >> x0=P.StartingPoint();
% >> d=inexactsubgradient(x,F,eps,0);
%
% Example 2; evaluate an inexact subgradient (absolute accuracy) of a
% Function F at some point x0
%
% >> P=pet(); param.L=1;param.mu=0.2; eps=.1;
% >> F=P.AddObjective('SmoothStronglyConvex',param);
% >> x0=P.StartingPoint();
% >> d=inexactsubgradient(x,F,eps,1);
%
d=Point('Point');
[g,~]=f.oracle(x);
if nargin < 4 % relative accuracy
    abs=0;
end

if abs==1 % absolute accuracy
    f.AddConstraint((g-d)^2-eps^2<=0);
else
    f.AddConstraint((g-d)^2-eps^2*g^2<=0);
end


end

