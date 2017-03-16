function [ g ] = subgradient( x,f )
% [ g ] = subgradient( x,f )
%
% This routine allows to evaluate a (sub)gradient. 
%
% Inputs:
%   - x is a point in the decision space (not a function value)
%   - f is a function
% Output:
%   - g is a subgradient of f at x
%
% Example; evaluate a subgradient of a Function F at some point x0
%
% >> P=pet(); param.L=1;param.mu=0.2; eps=.1;
% >> F=P.AddObjective('SmoothStronglyConvex',param);
% >> x0=P.StartingPoint();
% >> g=subgradient(x,F);
%
% Another possibility for evaluating (sub)gradients is to use the routine
% oracle featured by functions:
%
% >> P=pet(); param.L=1;param.mu=0.2; eps=.1;
% >> F=P.AddObjective('SmoothStronglyConvex',param);
% >> x0=P.StartingPoint();
% >> [g,f]=F.oracle(x);
%
% with f=F(x) the corresponding function value.
%
[g,~]=f.oracle(x);
end

