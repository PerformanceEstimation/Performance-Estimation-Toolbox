function A_ProximalPointMethod
% In this example, we use a proximal point method for solving
% a monotone inclusion problem
%   find x st   0 \in Ax
% where A is maximally monotone. We denote by JA the resolvents of A.
%
% We show how to compute the worst-case of ||x_{N}-x_{N-1}||^2 given that
% ||x_0-x_*||^2 <= 1 (where x_* is such that 0\in Ax_*). The result can be
% compared to
%
% [1] Guoyong Gu, and Junfeng Yang. "Optimal nonergodic sublinear
%     convergence rate of proximal point algorithm for maximal monotone
%     inclusion problems." (2019)


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the class of monotone inclusions
A = P.DeclareFunction('Monotone');

% (2) Set up the starting points
x0 = P.StartingPoint();
xs = A.OptimalPoint();


P.InitialCondition((x0-xs)^2<=1);  % Normalize the initial distance ||w0-ws||^2 <= 1

% (3) Algorithm
alpha = 2;
N     = 10;

x    = cell(N+1,1);
x{1} = x0;
for i = 1:N
    x{i+1} = proximal_step(x{i},A,alpha);
end


% (4) Set up the performance measure: ||x{i+1}-x{i}||^2
P.PerformanceMetric((x{N+1}-x{N})^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output: result should be: (1-1/N)^(N-1) / N, see [1]
[ double((x{N+1}-x{N})^2) (1-1/N)^(N-1) / N]
end