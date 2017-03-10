clear all; clc;
% In this example, we use a conditional gradient method for
% solving the constrained smooth convex minimization problem
%   min_x F(x)=f_1(x)+f_2(x); 
%   for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is L-smooth and convex and where f_2(x) is
% a convex indicator function of diameter at most D.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with a feasible point.


% (0) Initialize an empty PEP
P=pet();

% (1) Set up the objective function
paramf1.L=1;      % Smoothness parameter
paramf2.D=1;      % Diameter of the indicator function

f1=P.AddObjective('SmoothStronglyConvex',paramf1); 
f2=P.AddObjective('ConvexIndicator',paramf2);
F=f1+f2;

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		% x0 is some starting point
[xs,fs]=F.OptimalPoint();   % xs is an optimal point, and fs=F(xs)
[g0, f0]=f1.oracle(x0);
[s0, h0]=f2.oracle(x0);     % this makes x0 feasible (we impose the 
                            % existence of a subgradient s0 of f2 at x0).

x=x0;
N=10;
for i=1:N
    [dir,~]=f1.oracle(x);
    y=linearoptimization_step(dir,f2);
    lambda=2/(1+i);
    x=(1-lambda)*x+lambda*y;
end

% (4) Set up the performance measure
[g,f]=F.oracle(x);         % g=grad F(x), f=F(x)
P.PerformanceMetric(f-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% Should be better than the standard guarantee from Conditional Gradient:
% 2*paramf1.L*paramf2.D^2/(N+2)
%
% See Jaggi, Martin. "Revisiting Frank-Wolfe: Projection-free sparse 
%     convex optimization." In: Proceedings of the 30th International
%     Conference on Machine Learning (ICML-13), pp. 427–435 (2013)