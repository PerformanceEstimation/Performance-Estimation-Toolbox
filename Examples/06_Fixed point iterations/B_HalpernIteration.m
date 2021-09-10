function B_HalpernIteration
% In this example, we use the Halpern iteration for finding a fixed point
% to the non-expansive operator A :
%   find x such that  x = Ax
%
% [1] Lieder, Felix. "On the Convergence Rate of the Halpern-Iteration." 
%     (2017)


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
paramA.L=1;      % A is 1-Lipschitz (non-expansive)
A = P.DeclareFunction('Lipschitz',paramA);

% (2) Set up the starting point and initial condition
x0 = P.StartingPoint();		 % x0 is some starting point
xs = fixedpoint(A);
P.InitialCondition((x0-xs)^2<=1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 10;
lambda = @(k)(1/(k+2));
x=x0;
for i=1:N
    x = lambda(i-1) * x0 + (1-lambda(i-1)) * A.evaluate(x);
end
xN  = x;
AxN = A.evaluate(xN);
% (4) Set up the performance measure
P.PerformanceMetric((xN-AxN)^2); % Worst-case squared residual

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double((xN-AxN)^2)   % worst-case squared residual

% The result should be (2/(N+1))^2 as in [1]
end