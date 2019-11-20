clear all; clc;
% In this example, we use SAGA for solving the finite sum minimization
% problem
%   min_x {F(x)= 1/n(f1(x)+ ... + fn(x)) }
% for notational convenience we denote xs=argmin_x F(x).
% Functions f1,...,fn are assumed L-smooth and mu-strongly convex, and with
% proximal operators available.
%
% This code compute the exact rate for the Lyapunov function from the
% original SAGA paper:
% [1] Aaron Defazio. "A Simple Practical Accelerated Method for Finite
% Sums." (2014)
% (Theorem 5 of [1])

% (0) Initialize an empty PEP
P = pep();

% Problem setup:
L = 1;	% Smoothness
m = .1; % Strong convexity
n = 10; % Number of functions in the finite sum

gamma = sqrt((n-1)^2+4*n*L/m)/2/L/n-(1-1/n)/2/L; % stepsize
c     = 1/m/L;                                   % Lyapunov parameter
kappa = m*gamma/(1+m*gamma);                     

% (1) Set up the objective function
param.L   = L;
param.mu  = m;

% Define the objective function F:
f{1} = P.DeclareFunction('SmoothStronglyConvex',param);
F    = f{1}/n;             
for i = 2:n
    f{i} = P.DeclareFunction('SmoothStronglyConvex',param);
    F    = F + f{i}/n;
end
[xs, ~] = F.OptimalPoint('opt');    % this is x*

% (2) Initial values
gk = cell(n,1);
for i = 1:n
    gk{i} = P.StartingPoint();
end
x{1} = P.StartingPoint();           % this is xk

% (3) Initial value of the Lyapunov function
T0 =  (xs-x{1})^2;
for i = 1:n
    [gis, fis] = f{i}.oracle('opt');
    T0  = T0 + c/n * (gis-gk{i})^2;
end

P.InitialCondition( T0 <= 1);

% (4) Explicit computation of the expected value of the Lyapunov function 
% after one iteration (so: expectation over n possible scenarios: one for 
% each element fi in the function).

T1avg = 0;
for i = 1:n
    w{i} = x{1} + gamma * gk{i};
    for j = 1:n
        w{i} = w{i} - gamma/n * gk{j};
    end
    [x{i+1},gx] = proximal_step( w{i}, f{i},  gamma);
    T1 =  (xs-x{1+i})^2;
    for j = 1:n
        [gjs, fjs] = f{j}.oracle('opt');
        if i ~= j
            T1  = T1 + c/n * (gk{j}-gjs)^2;
        else
            T1  = T1 + c/n * (gjs-gx)^2;
        end
    end
    T1avg = T1avg + T1/n;
end
P.PerformanceMetric(T1avg);

% (5) Solve the PEP
P.solve()

% PESTO output vs Theorem 5 of [1]
[double(T1avg) (1-kappa)]













