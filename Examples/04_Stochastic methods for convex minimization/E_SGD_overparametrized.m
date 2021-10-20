function E_SGD_overparametrized
% In this example, we use stochastic gradient descent for solving the
% finite sum minimization problem
%   min_x {F(x)= 1/n [f1(x)+ ... + fn(x)]  }
% for notational convenience we denote xs=argmin_x F(x).
% Functions f1,...,fn are assumed L-smooth and mu-strongly convex.
%
% In addition, we assume a bounded variance at optimality:
%           E||fi'(x*)||^2 == 0, (i.e., fi'(x*) = 0)
% which is standard from the SGD literature.
%
% This code compute the exact worst-case guarantee for the distance to
% optimality. We will observe it does not depend on n for this particular
% setting, meaning that the guarantees are also valid for expectation
% minimization settings (i.e., when n goes to infinity).
% 

% Initialize an empty PEP
P = pep();

% Problem setup:
L       = 1;        % Smoothness
m       = .1;       % Strong convexity
n       = 5;        % Number of functions in the finite sum
v       = 0;        % variance: E[||fi'(x*)||^2] <= v^2
gamma   = 1/L;      % stepsize
R       = 2;        % initial ||x0-xs||^2
kappa   = L/m;
% (1) Set up the objective function
param.L   = L;
param.mu  = m;

f{1} = P.DeclareFunction('SmoothStronglyConvex',param);
F = f{1}/n;              % F is the objective function
for i = 2:n
    f{i} = P.DeclareFunction('SmoothStronglyConvex',param);
    F    = F + f{i}/n;
end

x{1} = P.StartingPoint();
[xs, ~] = F.OptimalPoint('opt');

% Bounded variance:
var = 0;
for i = 1:n
    [gis, fis] = f{i}.oracle('opt');
    var        = var + gis^2/n;
end
P.InitialCondition( var <= v^2);
P.InitialCondition( (x{1}-xs)^2 <= R^2 );

% One iteration

distavg = 0;
for i = 1:n
    x{i+1}  = x{1} - gamma * f{i}.gradient(x{1});
    distavg = distavg + (x{i+1}-xs)^2/n;
end
P.PerformanceMetric(distavg);

% Solve the PEP
P.solve()

% should be equal to
[double(distavg) 1/2*(1-1/kappa)^2*R^2+1/2*(1-1/kappa)*R*sqrt((1-1/kappa)^2*R^2+4*v^2/L^2)+v^2/L^2]
end