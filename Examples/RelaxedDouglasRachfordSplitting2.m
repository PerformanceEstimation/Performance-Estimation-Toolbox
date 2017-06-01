clear all; clc;
% In this example, we use a Douglas-Rachford splitting (DRS) 
% method for solving the composite convex minimization problem
%   min_x F(x)=f_1(x)+f_2(x) 
%   (for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) and f_2 may potentially both have strong convexity and/or
% smoothness characteristics, and are both convex.
%
% We show how to compute the worst-case value of ||wN-ws||^2 when wN is
% obtained by doing N steps of (relaxed) DRS starting with an initial 
% iterate w0 satisfying ||w0-ws||<=1, and ws is some point to which
% the iterates of DRS converge.
%
% The algorithm is as follows:
%       x_{k+1}=prox_{\lambda f}(w_k)
%       y_{k+1}=prox_{\lambda g}(2*x_{k+1}-w_k)
%       w_{k+1}=(1-\alpha) w_k + \alpha (2*y_{k+1}-2*x_{k+1}+w_k)
%
% Note that the point ws may be defined in the following way:
% let g1s and g2s be subgradients of respectively f_1 and f_2 at xs such 
% that g1s+g2s=0 (optimality conditions). Then, ws=xs+lambda*g2s (lambda is
% the step size used in DRS).

% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
paramf1.mu=.3;	% Strong convexity parameter of f1
paramf1.L=Inf;  % Smoothness parameter of f1
paramf2.mu=0;   % Strong convexity parameter of f2
paramf2.L=1;    % Smoothness parameter of f2
f1=P.DeclareFunction('SmoothStronglyConvex',paramf1);
f2=P.DeclareFunction('SmoothStronglyConvex',paramf2);
F=f1+f2; % F is the objective function

% (2) Set up the starting point and initial condition
w0=P.StartingPoint(); % x0 is some starting point
[xs,fs]=F.OptimalPoint('opt'); % xs is an optimal point, and fs=F(xs)

% note that we tag the point xs as 'opt' to be able to re-evaluate it
% easily (providing the oracle routine with this tag allows to recover
% previously evaluated points).

% the next step evaluates the oracle at the tagged point 'opt' (xs) for
% recovering the values of g1s and g2s; this allows to guarantee that
% g1s+g2s=0;
[g1s,~]=f1.oracle('opt');
[g2s,~]=f2.oracle('opt');
lambda=1; ws=xs+lambda*g2s;

% Add an initial condition ||w0-ws||^2<= 1
P.InitialCondition((w0-ws)^2-1<=0); 

% (3) Algorithm
N=1;            % number of iterations
gam=lambda;		% step size
alpha=1;      % relaxation parameter

w=w0;
for i=1:N
    x=proximal_step(w,f2,gam);
    y=proximal_step(2*x-w,f1,gam);
    w=(1-alpha)*w+alpha*(2*y-2*x+w);
end

% (4) Set up the performance measure
P.PerformanceMetric((w-ws)^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double((w-ws)^2)   % worst-case distance to fixed point ws

% The results should be compared to
%
% (1) case where f1 is strongly convex and smooth, and f2 is just convex
% (abs(1-alpha)+alpha*max((1-gam*paramf1.mu)/(1+paramf1.mu*gam),(gam*paramf1.L-1)/(1+gam*paramf1.L)))^(2*N) 
% within the assumptions of the following result:
%
% see (Theorem 2): Giselsson, Pontus, and Stephen Boyd. "Linear 
%                  convergence and metric selection in Douglas-Rachford
%                  splitting and ADMM."
%                  IEEE Transactions on Automatic Control (2016). 
%
% (2) case where f1 is strongly convex, and f2 smooth and convex; the bound
% obtained from PEP should be at most (under the assumption of the
% following theorem) :
%
% abs(1-2*alpha+alpha*(1/(gam*paramf1.mu)+gam*paramf2.L)/(1+1/(gam*paramf1.mu)+gam*paramf2.L))+alpha*(1/(gam*paramf1.mu)+gam*paramf2.L)/(1+1/(gam*paramf1.mu)+gam*paramf2.L)
%
% see (Theorem 5.6): Giselsson, Pontus. "Tight Global Linear Convergence
%                  Rate Bounds for Douglas-Rachford Splitting."
%                  Journal of Fixed Point Theory and Applications (2017). 




