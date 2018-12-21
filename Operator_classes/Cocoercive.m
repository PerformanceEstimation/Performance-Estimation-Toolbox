function cons=Cocoercive(pt1,pt2,beta)
%
% This routine implements the interpolation conditions for cocoercive
% operators
%
% To generate a 1-cocoercive operator 'h' from an instance of PEP called P:
%  >> P=pep();
%  >> param.beta = 1;
%  >> h=P.AddObjective('Cocoercive',param);
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

assert(beta>0,'Constants provided to the operator class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g))
    cons=((pt2.g-pt1.g)*(pt2.x-pt1.x)- beta *(pt2.g-pt1.g)^2>=0);
else
    cons=[];
end

end