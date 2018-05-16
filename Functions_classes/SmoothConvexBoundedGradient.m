function cons=SmoothConvexBoundedGradient(pt1,pt2,D,R,L)
%
% This routine implements the interpolation conditions for smooth convex
% functions with bounded subgradients. Three parameters may be provided: 
%       - a bound D on the maximum distance between any two subgradients 
%       (diameter-type bound: ||g_1-g_2||<= D for any g_1,g_2 being
%       subgradients of the function at (respectively) some x_1 and x_2)
%       (D is nonnegative, possibly infinite).
%       - a bound R on the maximum norm of any subgradient
%       (radius-type bound: ||g||<= R for any g being a
%       subgradient of the function at some x) (R is nonnegative, 
%       possibly infinite),
%       - smoothness constant L (L is nonnegative, possibly infinite).
%
% Note: not defining the value of L, D and/or R automatically corresponds 
%       to set L, D and/or R to infinity.
%
% To generate a smooth convex function 'f' with diameter-type bound
% D=infinite, radius-type bound R=1 and a smoothness constant L=1.5,
% from an instance of PEP called P:
%  >> P=pet();
%  >> param.D=Inf; param.R=1; param.L=1.5;
%  >> f=P.AddObjective('SmoothConvexBoundedGradient',param);
%
% For details about interpolation conditions, we refer to the following
% references:
%
% (1) Taylor, Adrien B., Julien M. Hendrickx, and François Glineur. 
%     "Smooth strongly convex interpolation and exact worst-case 
%     performance of first-order methods." 
%     Mathematical Programming 161.1-2 (2017): 307-345.
%
% (2) Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%     "Exact Worst-case Performance of First-order Methods for Composite
%     Convex Optimization."to appear in SIAM Journal on Optimization (2017)
%
%
assert(D>=0 & R>=0 & L>=0,'Constants provided to the functional class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    if L~=Inf
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+1/(2*L)*(pt1.g-pt2.g)^2)<=0);
    else
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x))<=0);
    end
    if D~=Inf
        cons=cons+((pt1.g-pt2.g)^2-D^2<=0);
    end
else
    if R~=Inf 
        cons=((pt1.g)^2-R^2<=0);
    else
        cons=[];
    end
end

end
