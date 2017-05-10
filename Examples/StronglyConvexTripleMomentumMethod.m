clear all; clc;
% In this example, we use the triple momentum method for solving the 
% L-smooth mu-strongly convex minimization problem
%   min_x F(x); 
%   for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.L=1;      % Smoothness parameter
param.mu=.1;    % Strong convexity parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();             % x0 is some starting point
[xs,fs]=F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
N=10;		% number of iterations

kappa=param.mu/param.L;
rho=(1-sqrt(kappa));

alpha=(1+rho)/param.L; beta=rho^2/(2-rho); gamma=rho^2/(1+rho)/(2-rho);
delta=rho^2/(1-rho^2);

xsi=cell(N+1,1);% we store the iterates in a cell for convenience
xsi{1}=x0;
xsi{2}=x0;
y=x0;
for i=2:N+1
    xsi{i+1}=(1+beta)*xsi{i}-beta*xsi{i-1}-alpha*F.gradient(y);
    y=(1+gamma)*xsi{i+1}-gamma*xsi{i};
    x=(1+delta)*xsi{i+1}-delta*xsi{i};
end

% (4) Set up the performance measure
fN=F.value(x);         % fN=F(xN)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% Should at most match the standard guarantees
% f(xN)-f(xs)<= (1-sqrt(kappa))^(2*N)*param.L/2/kappa*double((x0-xs)^2)
% ||xN-xs||^2<= (1-sqrt(kappa))^(2*N)/kappa*double((x0-xs)^2)
