function [x] = proximal_step(x0,func,gamma)
% [x] = proximal_step(x0,func,gamma)
%
% This routine performs a proximal step of step size gamma, starting from
% x0, and on function f. That is, it performs:
%       x=x0-gamma*g, where g is a (sub)gradient of f at x.
%       (implicit/proximal scheme).
%
% Input: - starting point x0
%        - function f on which the (sub)gradient will be evaluated
%        - step size gamma of the proximal step
%
% Output: x=x0-gamma*g, where g is a (sub)gradient of f at x.
%
g_imp=Point('Point');
x=x0-gamma*g_imp;
f=Point('Function value');
func.AddComponent(x,g_imp,f);

end

