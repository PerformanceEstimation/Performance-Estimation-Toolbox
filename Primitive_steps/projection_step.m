function x = projection_step(x0,ind,tag)
% [x] = projection_step(x0,ind)
%
% This routine performs a projection step of step size gamma, starting from
% x0, and on function f. That is, it performs:
%       x=Proj_ind(x0).
% Note that this is just equivalent to perform a proximal step on an
% indicator function.
%
% Input: - starting point x0
%        - indicator function of the set on which x0 should be projected
%        - optional tag
%
% Output: x=Proj_ind(x0).
%
g_ind=Point('Point');
x=x0-g_ind;
f=Point('Function value');
ind.AddComponent(x,g_ind,f);

if nargin > 2
    ind.AddComponent(x,g_ind,f,tag);
else
    ind.AddComponent(x,g_ind,f);
end
end

