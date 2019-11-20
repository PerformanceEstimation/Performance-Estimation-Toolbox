function cons=Smooth(pt1,pt2,L)
%
% This routine implements the interpolation conditions for non-convex
% smooth functions. One parameter may be provided: 
%       - smoothness constant L (L is nonnegative and possibly infinite).
%
% Note: not defining the value of L automatically corresponds to set L=Inf.
%
% To generate a non-convex smooth function 'f' with smoothness paramater
% L=1 from an instance of PEP called P:
%  >> P=pep();
%  >> param.L=1;
%  >> f=P.AddObjective('Smooth',param);
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
%     Convex Optimization." SIAM Journal on Optimization (2017)
%
%  *** There is a typo in [2] Theorem 6, in front of the scalair product
%  (should be a minus and not a plus); this is corrected below).
%
assert(L>=0,'Constants provided to the functional class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=((pt1.f-pt2.f-L/4*(pt1.x-pt2.x)^2-1/2*(pt1.g+pt2.g)*(pt1.x-pt2.x)+1/(4*L)*(pt1.g-pt2.g)^2)<=0);
else
    cons=[];
end

end