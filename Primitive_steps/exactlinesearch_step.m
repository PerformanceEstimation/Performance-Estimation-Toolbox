function [ x,gx,fx ] = exactlinesearch_step(x0,f,dirs,tag)
% [ x,gx,fx ] = exactlinesearch_step(x0,f,dirs)
%
% This routine MIMICS an exact line search.
%
% Input: - starting point x0
%        - function f on which the (sub)gradient will be evaluated
%        - directions dirs (cell) containing all directions required
%        to be orthogonal to the (sub)gradient of x. Note that (x-x0) is
%        automatically constrained to be orthogonal to the subgradient of f
%        at x.
%        - optional tag
%
% Output: - x such that all vectors in dirs are orthogonal to the gradient of
%         f at x.
%         - gx is a subgradient of f at x (the one satisfying optimality
%         conditions)
%         - fx is f(x)
%
% Example of its use in
%  "Examples/01_Methods for unconstrained convex minimization/F_GradientExactLineSearch.m"
%
x=Point('Point');
if nargin > 3
    [gx,fx]=f.oracle(x,tag);
else
    [gx,fx]=f.oracle(x);
end
nb_orth=max(size(dirs));
f.AddConstraint((x-x0)*gx==0,'Orthogonality condition (exact linesearch)');
if isa(dirs,'cell')
    for i=1:nb_orth
        f.AddConstraint(dirs{i}*gx==0,'Orthogonality condition (exact linesearch)');
    end
else
    f.AddConstraint(dirs*gx==0,'Orthogonality condition (exact linesearch)');
end
end

