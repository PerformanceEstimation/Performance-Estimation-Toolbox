function P_InexactFastGradientMethod
% In this example, we use an inexact fast gradient method for solving the 
% L-smooth convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x);
% where F(x) is L-smooth and convex. The inexactness model is in terms of
% relative inaccuracy:
%
%   (gi is the gradient of F at xi)
%       ||d-gi||<=eps*||gi||
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.mu = 0;      % strong convexity parameter
param.L  = 1;      % Smoothness parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0       = P.StartingPoint();        % x0 is some starting point
[xs, fs] = F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2 <= 1);  % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N = 5; % number of iterations

x    = cell(N+1,1); % store iterates in a cell
x{1} = x0;
y    = x0;
eps  = 0;
for i = 1:N
    d      = inexactsubgradient(y, F, eps);
    x{i+1} = y - 1/param.L * d;
    y      = x{i+1} + (i-1)/(i+2) * (x{i+1} - x{i});
end

% (4) Set up the performance measure
[g, f] = F.oracle(x{N+1});   
P.PerformanceMetric(f - fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
out=P.solve()

% (6) Evaluate the output
double(f - fs)   % worst-case objective function accuracy

% Result should be worse than 2/(N^2+5*N+6) in the case eps = 0
% (exact fast gradient)
%
% see Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%     "Exact Worst-case Performance of First-order Methods for Composite
%     Convex Optimization." SIAM Journal on Optimization (2017)

end



