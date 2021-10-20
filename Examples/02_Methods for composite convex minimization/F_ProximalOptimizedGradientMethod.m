function E_FastProximalGradientMethod
% In this example, we use the optimized gradient method (OGM) for
% solving the L-smooth convex minimization problem
%   min_x {F(x) = f1(x) + f2(x) }
% (for notational convenience we denote xs=argmin_x F(x)) with f1(x) being
% smooth and convex and f2 closed convex and proper with proximal operator
% available.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of proximal optimized gradient method (POGM)
% starting with an initial iterate satisfying ||x0-xs||<=1.
%
% Note that OGM is developped in the following two works:
%
% [1] Drori, Yoel, and Marc Teboulle.
%     "Performance of first-order methods for smooth convex minimization:
%     a novel approach." Mathematical Programming 145.1-2 (2014): 451-482.
%
% [2] Kim, Donghwan, and Jeffrey A. Fessler.
%     "Optimized first-order methods for smooth convex minimization."
%     Mathematical programming 159.1-2 (2016): 81-107.
% 
% whereas POGM was developed in
%
% [3]  Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%      "Exact Worst-case Performance of First-order Methods for Composite
%      Convex Optimization." SIAM Journal on Optimization (2017)
%


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
L = 1; param.L = L;      % Smoothness parameter

% F is the objective function:
f1 = P.DeclareFunction('SmoothStronglyConvex',param);
f2 = P.DeclareFunction('Convex',param);
F  = f1 + f2;

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
lambda = 1/param.L;  % stepsize
gamma  = 1;          % any nonzero initialization (first one is not used, see [3])
N      = 10;         % number of iterations

y        = cell(N+1,1);% store iterates in a cell
y{1}     = x0;
x        = x0;
z        = x0;
theta    = cell(N+1,1);
theta{1} = 1;
for i = 1:N
    y{i+1} = gradient_step(x, f1, lambda);
    if i<N
        theta{i+1}  = (1+sqrt(4*theta{i}^2+1))/2;
    else
        theta{i+1}  = (1+sqrt(8*theta{i}^2+1))/2;
    end
    z = y{i+1} + (theta{i}-1)/theta{i+1} * (y{i+1} - y{i}) ...
        + theta{i}/theta{i+1} * (y{i+1} - x) ...
        + (theta{i} - 1)/(L*gamma*theta{i+1}) * ( z - x);
    
    gamma = (2*theta{i}+theta{i+1}-1)/(L*theta{i+1});
    x     = proximal_step(z, f2, gamma);
end

% (4) Set up the performance measure
fN = F.value(x);              
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% The result should be larger than 1/2/theta{N+1}^2, while decreasing
% at the same speed, see [3].
end
