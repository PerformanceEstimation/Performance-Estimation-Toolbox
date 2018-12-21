function cons=LipschitzStronglyMonotone(pt1,pt2,L,mu)
%
% This routine implements quadratic constraints (no interpolation 
% conditions) operators for operators that are both Lipschitz and 
% strongly monotone.
%
%  **********************************************************************
%  NOTE: No interpolation conditions are known for this type of operators;
%  results are therefore NOT guaranteed to be tight.
%  **********************************************************************
%
% To generate a 1-Lipschitz .1-strongly monotone operator 'h' from 
% an instance of PEP called P:
%  >> P=pep();
%  >> param.L = 1; param.mu = .1;
%  >> h=P.AddObjective('LipschitzStronglyMonotone',param);
%
%  **********************************************************************
%  NOTE: PESTO was initially though for evaluating performances of
%  optimization algorithms. Operators are represented in the same way as
%  functions, but function values are not accessible.
%  **********************************************************************
%
% For details about interpolation conditions, we refer to the following:
%
% (1) E. K. Ryu, A. B. Taylor, C. Bergeling, and P. Giselsson, 
% "Operator Splitting Performance Estimation: Tight contraction factors 
%  and optimal parameter selection," arXiv:1812.00146, 2018.
%
%

assert(L>0 & mu>=0 & L>=mu,'Constants provided to the operator class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g))
    cons=((pt2.g-pt1.g)*(pt2.x-pt1.x)- mu *(pt2.x-pt1.x)^2>=0);
    cons= cons + ((pt2.g-pt1.g)^2- L^2 *(pt2.x-pt1.x)^2<=0);
else
    cons=[];
end

end