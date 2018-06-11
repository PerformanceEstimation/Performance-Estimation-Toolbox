function cons=SmoothStronglyConvex(pt1,pt2,mu,L)
%
% This routine implements the interpolation conditions for smooth strongly 
% convex functions. Two parameters may be provided: 
%       - smoothness constant L (L is nonnegative, possibly infinite).
%       - strong convexity constant mu (mu is nonnegative).
%
% Note: not defining the value of L automatically corresponds to set L to
%       infinity; not defining the value of mu automatically corresponds to
%       set mu to zero.
%
% To generate a smooth strongly convex function 'f' with smoothness 
% constant L=1.5 and strong convexity constant mu=0.1
% from an instance of PEP called P:
%  >> P=pep();
%  >> param.mu=0.1; param.L=1.5;
%  >> f=P.AddObjective('SmoothStronglyConvex',param);
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
assert(mu>=0 & L>=0 & L>=mu,'Constants provided to the functional class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    if L~=Inf
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+...
            1/(2*(1-mu/L))*(1/L*(pt1.g-pt2.g)*(pt1.g-pt2.g)+...
            mu*(pt1.x-pt2.x)*(pt1.x-pt2.x)-...
            2*mu/L*(pt1.x-pt2.x)*(pt1.g-pt2.g)))<=0);
    else
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+mu/2*(pt1.x-pt2.x)^2)<=0);
    end    
else
    cons=[];
end
end