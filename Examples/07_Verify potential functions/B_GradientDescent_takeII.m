function B_GradientDescent_takeII
% In this example, we use a fixed-step gradient method for
% solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to verify the inequality (see [1] below)
%  
%      phi_{k+1} <=  phi_k
%
% where phi_{k} = (2k+1)L (F(x_{k+1})-F(xs)) + k(k+2) ||F'(x_{k+1})||^2 
%                           + L/2 * ||x_{k+1}-xs||^2 
% 
% when x_{k+1} = x_k - 1/L F'(x_k).
%
% This potential is proved in [1].
%
% [1] Adrien Taylor, and Francis Bach. "Stochastic first-order
%     methods: non-asymptotic and computer-aided analyses via
%     potential functions." (2019)

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.mu = 0;     % Strong convexity parameter
L        = 1;
param.L  = L;     % Smoothness parameter

% F is the objective function
F = P.DeclareFunction('SmoothStronglyConvex',param); 

% (2) Set up the starting point 
xk      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)


% (3) Algorithm
h = 1/param.L;		% step size
k = 10;             % pick iteration counter

[  gk, fk ]   = F.oracle(xk);
xkp1          = xk - h * gk; % x_{k+1}
[gkp1, fkp1 ] = F.oracle(xkp1);

% (4) Set up the performance measure:
% find the largest possible value for phi_{k+1}-phi_k, which should be <= 0

phi=@(k,f,df,x)( L*(2*k+1)*f+k*(k+2)*df^2+L^2*x^2);
objective = phi(k+1,fkp1-fs,gkp1,(xkp1-xs))-phi(k,fk-fs,gk,(xk-xs));
P.PerformanceMetric(objective);    


% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(objective)   % worst-case objective function accuracy

% The result should be <= 0 (see [1])
end