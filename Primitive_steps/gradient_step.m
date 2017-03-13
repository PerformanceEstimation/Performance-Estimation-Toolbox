function [x] = gradient_step(x0,f,gamma)
% This routine performs a gradient step of step size gamma, starting from
% x0, and on function f. That is, it performs:
%       x=x0-gamma*g0, where g0 is a (sub)gradient of f at x0.
%
% Input: - starting point x0
%        - function f on which the (sub)gradient will be evaluated
%        - step size gamma of the gradient step
%
% Output: x=x0-gamma*g0, where g0 is a (sub)gradient of f at x0.
%
[g,~]=f.oracle(x0);
x=x0-gamma*g;

end

