function [ x ] = exactlinesearch_step(x0,f,dirs)
% [ x ] = exactlinesearch_step(x0,f,dirs)
%
% This routine MIMICS an exact line search.
%
% Input: - starting point x0
%        - function f on which the (sub)gradient will be evaluated
%        - directions dirs (cell) containing all the directions required to
%        be orthogonal to the (sub)gradient of x.
%
% Output: x such that all vectors in dirs are orthogonal to the gradient of
%         f at x.
%
% Example of its use in Examples/GradientExactLineSearch.m.
%
x=Point('Point');
[g,~]=f.oracle(x);
nb_orth=max(size(dirs));
f.AddConstraint((x-x0)*g<=0);
if isa(dirs,'cell')
    for i=1:nb_orth
        f.AddConstraint(dirs{i}*g<=0);
    end
else
    f.AddConstraint(dirs*g<=0);
end
end

