function [x] = proximal_step(x0,func,gamma)
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
%
% Output: x=x0-gamma*g, where g is a (sub)gradient of func at x.
%
g_imp=Point('Point');
x=x0-gamma*g_imp;
f=Point('Function value');
func.AddComponent(x,g_imp,f);

end

