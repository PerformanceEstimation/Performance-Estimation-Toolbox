function [x,gx,fx] = proximal_step(x0,func,gamma,tag)
% [x] = proximal_step(x0,func,gamma)
%
% This routine performs a proximal step of step size gamma, starting from
% x0, and on function func. That is, it performs:
%       x=x0-gamma*g, where g is a (sub)gradient of func at x.
%       (implicit/proximal scheme).
%
% Input: - starting point x0
%        - function func on which the (sub)gradient will be evaluated
%        - step size gamma of the proximal step
%        - optional tag
%
% Output: x=x0-gamma*g, where g is a (sub)gradient of func at x.
%
gx=Point('Point');
x=x0-gamma*gx;
fx=Point('Function value');

if nargin > 3
    func.AddComponent(x,gx,fx,tag);
else
    func.AddComponent(x,gx,fx);
end
end

