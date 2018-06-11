function cons=StronglyConvexBoundedDomain(pt1,pt2,D,R,mu)
%
% This routine implements the interpolation conditions for strongly convex
% functions with bounded domains. Three parameters may be provided: 
%       - a diameter D (D is nonnegative, and possibly infinite),
%       - a radius R (R is nonnegative, and possibly infinite),
%       - a strong convexity constant mu (mu is nonnegative).
%
% Note: not defining the value ofD and/or R automatically corresponds 
%       to set D and/or R to infinite. Not defining mu corresponds to set
%       it to zero.
%
% To generate a smooth convex function 'f' with diameter-type bound
% D=infinite, radius-type bound R=1 and strong convexity mu=1.5,
% from an instance of PEP called P:
%  >> P=pep();
%  >> param.D=Inf; param.R=1; param.mu=1.5;
%  >> f=P.AddObjective('StronglyConvexBoundedDomain',param);
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
    cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+mu/2*(pt1.x-pt2.x)^2)<=0);
    if D~=Inf && pt1.num<pt2.num
        cons=cons+((pt1.x-pt2.x)^2-D^2<=0);
    end
else
    if R~=Inf
        cons=((pt1.x)^2-R^2<=0);
    else
        cons=[];
    end
end

end