function [ d ] = inexactsubgradient( x,f,eps,abs )


d=Point('Point');
[g,~]=f.oracle(x);
if nargin < 4 % relative accuracy
    abs=0;
end

if abs==1 % absolute accuracy
    f.AddConstraint((g-d)^2-eps^2<=0);
else
    f.AddConstraint((g-d)^2-eps^2*g^2<=0);
end

% Another possibility:

% err=Point('Point');
% [g,~]=f.oracle(x);
% d=g+err;
% if nargin < 4 
%     abs=0;%=0: relative accuracy
% end
% 
% if abs==1
%     f.AddConstraint(err^2-eps^2<=0);
% else
%     f.AddConstraint(err^2-eps^2*g^2<=0);
% end

end

