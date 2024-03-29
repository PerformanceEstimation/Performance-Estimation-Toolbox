function B_ProximalPointMethod
% In this example, we use a proximal point method for solving the 
% non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
%
% We show how to obtain a low-dimensional problem satisfying 
%   F(x_N) - F(x_*) >= ||x_0-x_*||^2 / 4N
% for some values of N
%
% Alternative interpretations:
% (1) the following code compute the solution to the problem
%       max_{F,x0,...,xN,xs} (F(xN)-F(xs))/||x0-xs||^2 
%            s.t. x1,...,xN are generated via the proximal point method,
%                 F is closed, proper, and convex.
%     where the optimization variables are the iterates and the convex
%     function F.
%
% (2) the following code compute the smallest possible value of 
%     C(N, step sizes) such that the inequality
%       F(xN)-F(xs)  <= C(N, step sizes) * ||x0-xs||^2
%     is valid for any closed, proper and convex F and any sequence of
%     iterates x1,...,xN generated by the proximal point method on F.
%


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
F = P.DeclareFunction('Convex');      % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();        % x0 is some starting point
[xs,fs] = F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2  == 1);% Initial condition ||x0-xs||^2 == 1

% (3) Algorithm
N   = 10;		% number of iterations
gam = @(k)(1);	% step size (possibly a function of k)

x = x0;
for i = 1:N
    x = proximal_step(x, F, gam(i));
end
xN = x;

% (4) Set up the performance measure
fN = F.value(xN);
P.AddConstraint(fN-fs >= (x0-xs)^2 / 4/N);

% (5) Solve the PEP
P.TraceHeuristic(1); % tries to find a low-dimensional worst-case instance
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% The result should be (and is) 1/(4*\sum_{i=1}^N gam(i))
% see Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%     "Exact Worst-case Performance of First-order Methods for Composite
%     Convex Optimization.", SIAM Journal on Optimization (2017)
double(x0) % is this small-dimensional?
% this should be a scalar value (i.e., there exists a 1-D example)
end