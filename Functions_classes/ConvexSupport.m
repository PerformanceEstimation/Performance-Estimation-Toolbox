function cons=ConvexSupport(pt1,pt2,D,R)
%
% This routine implements the interpolation conditions for convex
% support functions. Two parameters may be provided: 
%       - a bound D on the maximum distance between any two subgradients 
%       (diameter-type bound: ||g_1-g_2||<= D for any g_1,g_2 being
%       subgradients of the function at (respectively) some x_1 and x_2)
%       (D is nonnegative, and possibly infinite).
%       - a bound R on the maximum norm of any subgradient
%       (radius-type bound: ||g||<= R for any g being a
%       subgradient of the function at some x)
%       (R is nonnegative, and possibly infinite).
%
% Note: not defining the value of D and/or R automatically corresponds to 
%       set D and/or R to Inf.
%
% To generate a convex support function 'h' with diameter-type bound
% D=infinite and a radius-type bound R=1 from an instance of PEP called P:
%  >> P=pep();
%  >> param.D=Inf; param.R=1;
%  >> h=P.AddObjective('ConvexSupport',param);
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
assert(D>=0 & R>=0,'Constants provided to the functional class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=(pt1.x*(pt2.g-pt1.g)<=0);
    if D~=Inf
        cons=cons+((pt1.g-pt2.g)^2-D^2<=0);
    end
else
    cons=(pt1.x*pt1.g-pt1.f==0);
    if R~=Inf
        cons=cons+((pt1.g)^2-R^2<=0);
    end
end
end
