function [x, sx, hx] = mirror(gx0, sx0, mirror_map, gamma, tag)
% [x,gx,hx] = mirror(gx0, sx0, func, gamma, tag)
%
% This routine performs a mirror step of step size gamma.
% That is, denoting by h(.) the mirror map, it performs
%  h'(x) = h'(x0) - gamma*gx0, where h'(x) and h'(x0) are respectively
%                             gradients of the mirror map at x and x0.
%
% NOTE: it assumes the mirror map is differentiable.
%
%
% Input: - starting gradient sx0 (e.g., gradient at x0 of 'mirror_map'),
%        - step gx0 (e.g., gradient at x0 of the function to be minimized),
%        - mirror_map on which the (sub)gradient will be evaluated,
%        - step size gamma,
%        - optional tag.
%
% Output: - x:  mirror point,
%         - sx: subgradient of the mirror_map at x that was used in the procedure,
%         - hx: value of the mirror_map evaluated at x.
%

hx = Point('Function value');
x  = Point('Point');

sx = sx0 - gamma * gx0;
if nargin > 4
    mirror_map.AddComponent(x,sx,hx,tag);
else
    mirror_map.AddComponent(x,sx,hx);
end
end

