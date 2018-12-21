
function cons=Monotone(pt1,pt2)
%
% This routine implements the interpolation conditions for monotone
% operators
%
% To generate a monotone operator 'h' from an instance of PEP called P:
%  >> P=pep();
%  >> h=P.AddObjective('Monotone');
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

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g))
    cons=((pt2.g-pt1.g)*(pt2.x-pt1.x)>=0);
else
    cons=[];
end

end