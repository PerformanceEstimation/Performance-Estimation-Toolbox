function x = fixedpoint(A)
% x = fixedpoint(A)
%
% This routine allows to find a fixed point of an operator A.
%
% Inputs: A is an operator or a function
%
% Output: 
%       - if A is an operator, x such that x = A x (fixed point of A)
%       - if A is a function, x such that x = A.gradient(x)
%
x=Point('Point');
f=Point('Function value');
A.AddComponent(x,x,f);


end

