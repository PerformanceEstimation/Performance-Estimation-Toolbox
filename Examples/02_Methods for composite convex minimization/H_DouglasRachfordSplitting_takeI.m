clear all; clc;
% In this example, we use a Douglas-Rachford splitting (DRS) 
% method for solving the composite convex minimization problem
%   min_x { F(x) = f_1(x) + f_2(x) }
%   (for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is L-smooth and mu-strongly convex, and f_2 is convex,
% closed and proper. Both proximal operators are assumed to be available.
%
% We show how to compute a contraction factor for the iterates of DRS 
% (i.e., how do the iterates contract when the algorithm is started from
% two different initial points).
%
%
% Our notations for the DRS algorithm are as follows:
%       x_k     = prox_{\alpha f2}(w_k)
%       y_k     = prox_{\alpha f1}(2*x_k-w_k)
%       w_{k+1} = w_k +\theta (y_k - x_k)
%
% and our goal is to compute the smallest contraction factor rho such that
% when the algorithm is started from two different points w_0 and w_0', we
% have ||w_1 - w_1'||^2 <= rho^2 ||w_0 - w_0'||^2.
%
% Details on the SDP formulations can be found in 
% [1] Ernest K. Ryu, Adrien B. Taylor, Carolina Bergeling,
%     and Pontus Giselsson. "Operator splitting performance estimation:
%     Tight contraction factors and optimal parameter selection." (2018)
%
% When theta = 1, the bound can be compared with that of 
% [2] Giselsson, Pontus, and Stephen Boyd. "Linear convergence and
%     metric selection in Douglas-Rachford splitting and ADMM."
%     IEEE Transactions on Automatic Control (2016). 

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
paramf1.mu = .1;        % Strong convexity parameter
paramf1.L  = 1;         % Smoothness parameter

f1 = P.DeclareFunction('SmoothStronglyConvex',paramf1);
f2 = P.DeclareFunction('Convex');
F  = f1 + f2; % F is the objective function

% (2) Set up the starting points and initial condition
w0   = P.StartingPoint(); % w0 is some starting point
w0p  = P.StartingPoint(); % w0' is some starting point

P.InitialCondition((w0-w0p)^2<=1); % Initial condition ||w0-w0'||^2<= 1

% (3) Algorithm
N       = 3;    % number of iterations
alpha   = 3;    % step size
theta   = 1;    % overrelaxation

% compute trajectory starting from w0
w = w0;
for i = 1:N
    x = proximal_step(w,f2,alpha);
    y = proximal_step(2*x-w,f1,alpha);
    w = w + theta * (y-x);
end
% compute trajectory starting from w0'
wp = w0p;
for i = 1:N
    xp = proximal_step(wp,f2,alpha);
    yp = proximal_step(2*xp-wp,f1,alpha);
    wp = wp + theta * (yp-xp);
end
% (4) Set up the performance measure (how do the trajectories contract?)
P.PerformanceMetric((w-wp)^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double((w-wp)^2)   % worst-case distance to fixed point ws

% When theta = 1, the result should be (and is) (see [2])
% max(1/(1+paramf1.mu*alpha),alpha*paramf1.L/(1+alpha*paramf1.L))^(2*N)






