function cons=ConvexIndicator(pt1,pt2,D,R)
%
% This routine implements the interpolation conditions for convex
% indicator functions with bounded domains. Two parameters may be provided: 
%       - a diameter D (D is nonnegative, and possibly infinite),
%       - a radius R (R is nonnegative, and possibly infinite).
%
% Note: not defining the value of D and/or R automatically corresponds to 
%       set D and/or R to infinite.
%
% To generate a convex indicator function 'h' with diameter D=infinite and
% a radius R=1 from an instance of PEP called P:
%  >> P=pet();
%  >> param.D=Inf; param.R=1;
%  >> h=P.AddObjective('ConvexIndicator',param);
%
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
if ~(pt1.x.isEqual(pt2.x))
    cons=(pt1.g*(pt2.x-pt1.x)<=0);
    if D~=Inf
        cons=cons+((pt1.x-pt2.x)^2-D^2<=0);
    end
else
    cons=(pt1.f==0);
    if R~=Inf
        cons=cons+((pt1.x)^2-R^2<=0);
    end
end

end