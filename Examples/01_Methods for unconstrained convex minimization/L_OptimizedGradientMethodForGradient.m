function L_OptimizedGradientMethodForGradient
% In this example, we use the optimized gradient method for gradient norm
% (OGM-G) for solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of ||F'(xN)||^2 when xN is
% obtained by doing N steps of OGM-G starting with an initial
% iterate satisfying F(x0)-F(x*)<=1.
%
% Note that OGMG is developped in the following work:
% [1] Donghwan Kim, and Jeffrey Fessler. "Optimizing the Efficiency of 
%     First-order Methods for Decreasing the Gradient of Smooth 
%     Convex Functions." (2018)
%


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L = 1;      % Smoothness parameter

% F is the objective function
F = P.DeclareFunction('SmoothStronglyConvex',param);

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();    % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 	% xs is an optimal point, and fs=F(xs)
[g0,f0] = F.oracle(x0);
P.InitialCondition(f0-fs <= 1); % Initial condition F(x0)-F(x*)<=1.

% (3) Algorithm
N = 10;		% number of iterations

x = cell(N+1,1);% store iterates in a cell
y = cell(N,1);  % store iterates in a cell
g = cell(N+1,1);% store gradients in a cell
f = cell(N+1,1);% store function values in a cell
x{1} = x0; g{1} = g0; y{1} = x0;

theta    = cell(N+1,1);
theta{1} = 1;
for i=1:N
    if i<N
        theta{i+1} = (1+sqrt(4*theta{i}^2+1))/2;
    else
        theta{i+1} = (1+sqrt(8*theta{i}^2+1))/2;
    end
end
theta_inverse = @(i)theta{N+1-i};

for i = 1:N
    y{i+1} = x{i} - 1/param.L * g{i};
    cc     = (2*theta_inverse(i)-1)/(2*theta_inverse(i-1)-1);
    x{i+1} = y{i+1} + (theta_inverse(i-1)-1)/theta_inverse(i-1)*cc*(y{i+1}-y{i}) ...
        + cc * (y{i+1}-x{i});
    g{i+1} = F.gradient(x{i+1});
end

% (4) Set up the performance measure
P.PerformanceMetric(g{N+1}^2); % Worst-case performance of ||F'(xN)||^2

% (5) Solve the PEP
P.solve();

% (6) Evaluate the output
double(g{N+1}^2)   % worst-case value of ||F'(xN)||^2

% The result should be 2/theta{N+1}^2 (see [1])
end
