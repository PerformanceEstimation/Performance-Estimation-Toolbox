clear all; clc;
% In this example, we use a projected gradient method for
% solving the constrained smooth strongly convex minimization problem
%   min_x F(x)=f_1(x)+f_2(x); 
%   for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is L-smooth and mu-strongly convex and where f_2(x) is
% a convex indicator function.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying F(x0)-F(xs)<=1.

% (0) Initialize an empty PEP
P=pet();

% (1) Set up the objective function
paramf1.mu=.1;	% Strong convexity parameter
paramf1.L=1;    % Smoothness parameter
f1=P.AddObjective('SmoothStronglyConvex',paramf1);
f2=P.AddObjective('ConvexIndicator');
F=f1+f2; % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 	 % xs is an optimal point, and fs=F(xs)
[g0,f0]=F.oracle(x0);
P.InitialCondition(f0-fs<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gam=1/paramf1.L;		% step size
N=1;		% number of iterations

x=x0;
for i=1:N
    xint=gradient_step(x,f1,gam);
    x=projection_step(xint,f2);
end
xN=x;
fN=F.value(xN);

% (4) Set up the performance measure
P.PerformanceMetric(fN-fs);

% (5) Solve the PEP
P.solve()

% Result should be (and is) max((1-paramf1.mu*gam)^2,(1-paramf1.L*gam)^2)
