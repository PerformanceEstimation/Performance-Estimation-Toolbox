function [x, sx, hx, gx, fx] = mirror_prox(sx0, mirror_map, min_function, gamma,tag)
% [x,gx,hx] = mirror_prox(sx0, mirror_map, min_function, gamma, tag)
%
% This routine performs a proximal mirror step of step size gamma.

% That is, denoting by h(.) the mirror map, and f(.) the function to be 
% minimized, it performs
%  h'(x) = h'(x0) - gamma*f'(x), where h'(x) and h'(x0) are respectively
%                             (sub)gradients of the mirror map at x and x0,
%                             and f'(x) is a subgradient of f at x.
%
% Input: - starting gradient sx0 (e.g., gradient at x0 of 'mirror_map'),
%        - mirror_map on which the (sub)gradient will be evaluated,
%        - min_function which we aim to minimize,
%        - step size gamma,
%        - optional tag.
%
% Output: - x:  mirror point,
%         - sx: subgradient of the mirror_map at x that was used in the procedure,
%         - hx: value of the mirror_map evaluated at x,
%         - sx: subgradient of the min_function at x that was used in the procedure,
%         - hx: value of the min_function evaluated at x.

hx = Point('Function value');
x  = Point('Point');
fx = Point('Function value');
gx  = Point('Point');

sx = sx0 - gamma * gx;
if nargin > 4
    mirror_map.AddComponent(x,sx,hx,tag);
    min_function.AddComponent(x,gx,fx,tag);
else
    mirror_map.AddComponent(x,sx,hx);
    min_function.AddComponent(x,gx,fx);
end
end

