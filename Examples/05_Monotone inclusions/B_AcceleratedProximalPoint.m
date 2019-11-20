clear all; clc;
% In this example, we use an accelerated proximal point method for solving
% a monotone inclusion problem
%   find x st   0 \in Ax
% where A is maximally monotone. We denote by JA the resolvents of A.
%
% We show how to compute the worst-case of ||x_{N}-x_{N-1}||^2 given that
% ||x_0-x_*||^2 <= 1 (where x_* is such that 0\in Ax_*). The result can be
% compared to
%
% [1] Donghwan Kim. "Accelerated Proximal Point Method  and Forward Method
%     for Monotone Inclusions." (2019)


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
y    = cell(N+1,1);
x{1} = x0;
y{1} = x0;
y{2} = x0;

for i = 0:N-1
    x{i+2} = proximal_step(y{i+2},A,alpha);
    y{i+3} = x{i+2} + i/(i+2) * (x{i+2}-x{i+1}) - i/(i+2) * (x{i+1} - y{i+1});
end


% (4) Set up the performance measure:
P.PerformanceMetric((x{N+1}-y{N+1})^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output: result should be: 1/N^2, see [1]
[ double((x{N+1}-y{N+1})^2) 1/N^2]

