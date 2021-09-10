function C_ProximalGradientMethod
% In this example, we use a proximal gradient method for
% solving the composite convex minimization problem
%   min_x {F(x) = f_1(x)+f_2(x) } 
%   for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is L-smooth and mu-strongly convex and where f_2(x) is
% closed convex and proper.
%
% We show how to compute the worst-case value of ||xN-xs||^2 when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
paramf1.mu = .1;        % Strong convexity parameter
paramf1.L  = 1;         % Smoothness parameter

f1 = P.DeclareFunction('SmoothStronglyConvex',paramf1);
f2 = P.DeclareFunction('Convex');
F  = f1 + f2; % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();            % x0 is some starting point
[xs,fs] = F.OptimalPoint();             % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1);    % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
gamma = 1/paramf1.L;    % stepsize
N     = 1;                    % number of iterations

x = x0;
for i = 1:N
    xint = gradient_step(x, f1, gamma);
    x    = proximal_step(xint, f2, gamma);
end
xN = x;

% (4) Set up the performance measure
P.PerformanceMetric((xN-xs)^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double((xN-xs)^2)   % worst-case objective function accuracy

% Result should be max((1-paramf1.mu*gamma)^2,(1-paramf1.L*gamma)^2)

end


