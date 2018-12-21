clear all; clc;
% In this example, we use the optimized gradient method for gradient norm
% (OGM-G) for solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of ||F'(xN)||^2 when xN is
% obtained by doing N steps of OGM-G starting with an initial
% iterate satisfying F(x0)-F(x*)<=1.
%
% Note that OGMG is developped in the following work:
% (1) Kim, D., & Fessler, J. A. (2018). "Optimizing the Efficiency of 
%     First-order Methods for Decreasing the Gradient of Smooth 
%     Convex Functions." preprint arXiv:1803.06600.
%


% (0) Initialize an empty PEP
P=pep();
L = 1;
% (1) Set up the objective function
param.L=L;      % Smoothness parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();           % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		% xs is an optimal point, and fs=F(xs)
[g0,f0]=F.oracle(x0);
P.InitialCondition(f0-fs<=1);   % Add an initial condition F(x0)-F(x*)<=1.

% (3) Algorithm
gam=1/param.L;		% step size
N=5;		% number of iterations

x=cell(N+1,1);%we store all the x's in a cell (for convenience)
x{1}=x0;
g{1}=g0;
y{1}=x0;
theta(1)=1;
for i=1:N
    if i<N
        theta(i+1)=(1+sqrt(4*theta(i)^2+1))/2;
    else
        theta(i+1)=(1+sqrt(4*theta(i)^2+1))/2;
    end
end
th_ti = @(i)theta(N+1-i);

for i = 1:N
    y{i+1} = x{i} - 1/L * g{i};
    cc       = (2*th_ti(i)-1)/(2*th_ti(i-1)-1);
    x{i+1} = y{i+1} + (th_ti(i-1)-1)/th_ti(i-1)*cc*(y{i+1}-y{i}) + cc * (y{i+1}-x{i});
    g{i+1}   = F.gradient(x{i+1});
end

% (4) Set up the performance measure
% gN=F.gradient(x{N+1});                % g=grad F(x), f=F(x)
obj = (g{N+1})^2;
P.PerformanceMetric(obj); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve();

% (6) Evaluate the output
double(obj)   % worst-case objective function accuracy

% The result should be 2/theta(end)^2
% see: Kim, D., & Fessler, J. A. (2018). "Optimizing the Efficiency of 
%      First-order Methods for Decreasing the Gradient of Smooth 
%      Convex Functions." preprint arXiv:1803.06600.