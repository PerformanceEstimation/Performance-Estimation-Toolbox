clear all; clc;
% In this example, we use a Douglas-Rachford splitting (DRS) 
% method for solving the composite convex minimization problem
%   min_x { F(x) = f_1(x) + f_2(x) }
%   (for notational convenience we denote xs=argmin_x F(x);
% where f_2 is smooth and convex, and f_1 is closed, proper and convex.
% Both proximal operators are assumed to be available.
%
% Our notations for the DRS algorithm is as follows:
%       x_k     = prox_{\alpha f2}(w_k)
%       y_k     = prox_{\alpha f1}(2*x_k-w_k)
%       w_{k+1} = w_k +\theta (y_k - x_k)
%
% We show how to compute the worst-case value of F(yN)-F(xs) when yN is
% obtained by doing N steps of (relaxed) DRS starting with an initial 
% iterate w0 satisfying ||x0-xs||<=1.
% It is known that xk and yk converge to xs, but not wk, and hence 
% we require the initial condition on x0 (arbitrary choice; partially
% justified by the fact we choose f2 to be the smooth function).
% Note that yN is feasible as it has a finite value for f1 
% (output of the proximal opertor on f1) and as f2 is smooth.

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
paramf1.mu = 0;     % Strong convexity parameter of f1
paramf1.L  = Inf;   % Smoothness parameter of f1
paramf2.mu = 0;     % Strong convexity parameter of f2
paramf2.L  = 1;     % Smoothness parameter of f2
f1 = P.DeclareFunction('SmoothStronglyConvex',paramf1);
f2 = P.DeclareFunction('SmoothStronglyConvex',paramf2);
F  = f1 + f2; % F is the objective function

% (2) Set up the starting point
w0      = P.StartingPoint();    % w0 is some starting point
[xs,Fs] = F.OptimalPoint();     % xs is an optimal point, and fs=F(xs)

% (3) Algorithm
N     = 5;      % number of iterations
alpha = 1;		% step size
theta = 1;      % overrelaxation parameter

w     = cell(N+1,1);
x     = cell(N,1);
y     = cell(N,1);
fyval = cell(N,1);
w{1}  = w0;
for i = 1:N
    x{i}                = proximal_step(w{i},f2,alpha);
    [y{i},~,fyval{i}]   = proximal_step(2*x{i}-w{i},f1,alpha);
    w{i+1}              = w{i} + theta * (y{i}-x{i});
end

P.InitialCondition( (x{1}-xs)^2 <= 1); % Initial condition ||x0-xs||^2 <= 1

% (4) Set up the performance measure
F_final = f2.value(y{N})+fyval{N};
P.PerformanceMetric(F_final-Fs);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(F_final-Fs)   

% For convenience, here are a few outputs from the toolbox (illustrating
% the expected 1/k behavior of function values)
% N_list       = 1:10;
% pesto_output = [1/4 0.1273 0.0838 0.0627 0.0501 0.0417...
%     0.0357 0.0313 0.0278 0.0250];
% close all; loglog(N_list,pesto_output,'-r'); 






