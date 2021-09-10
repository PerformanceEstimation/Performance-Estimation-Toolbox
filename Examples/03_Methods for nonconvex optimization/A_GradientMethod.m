function A_GradientMethod
% In this example, we use a fixed-step gradient method for
% solving the L-smooth minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst value of min_{0<=i<=N} ||F'(xi)||^2
% when xN is obtained by doing N steps of the gradient method and when
% f(x0) - f(xN) <= 1.


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L  = 1;     % Smoothness parameter

% F is the objective function
F = P.DeclareFunction('Smooth',param);


% (3) Algorithm
h = 1/param.L;		% step size
N = 2;             % number of iterations

x       = cell(N+1,1);
g       = cell(N+1,1);
f       = cell(N+1,1);

x{1}        = P.StartingPoint();		 % x0 is some starting point
[g{1},f{1}] = F.oracle(x{1});
[xs, Fs]= F.OptimalPoint();
for i = 1:N
    x{i+1} = x{i} - h * g{i};
    P.PerformanceMetric(g{i}^2);
    [g{i+1},f{i+1}] = F.oracle(x{i+1});
end

% (4) Set up the performance measure
P.InitialCondition( f{1} - f{N+1}<= 1);
P.PerformanceMetric(g{N+1}^2);     

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(g{N+1}^2)   % worst-case objective function accuracy
end
