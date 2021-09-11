function P_ImprovedInteriorGradientAlgorithm_noConstraint
% In this example, we use the improved interior gradient algorithm (IGA) 
% method for solving the composite convex minimization problem
%   min_x {F(x) = f1(x) + f2(x)}  
%   (for notational convenience we denote xs=argmin_x F(x);
% where f1 is L-smooth and convex and f2 is a closed convex indicator.
% We use a kernel h that is assumed to be smooth strongly convex (see [1])
%
% THIS CODE IS A SIMPLIFICATION OF
%   "N_ImprovedInteriorGradientAlgorithm.m"
% FOR WHEN f2 = 0 (i.e., no constraint). The computations are therefore 
% slightly lighter.
%
% The method originates from:
% [1] Alfred Auslender, and Marc Teboulle. "Interior gradient and proximal
%     methods for convex and conic optimization."
%     SIAM Journal on Optimization (2006).
%

% (0) Initialize an empty PEP
P = pep();


% (1) Set up the objective function
paramf1.L = 1;      % Smoothness parameter of function f1
paramh.L  = Inf;    % Smoothness of the kernel h
paramh.mu = 1;      % strong convexity parameter of the kernel h

f1 = P.DeclareFunction('SmoothStronglyConvex',paramf1); 
h  = P.DeclareFunction('SmoothStronglyConvex',paramh); 

F  = f1;
% (2) Set up the starting point and initial condition
x0        = P.StartingPoint();             % x0 is some starting point
[xs,fs]   = F.OptimalPoint();              % xs is an optimal point, and fs=F(xs)
[sxs,hxs] = h.oracle(xs);                  
[sx0,hx0] = h.oracle(x0);                  
[gx0,fx0] = f1.oracle(x0);

% (3) Algorithm
L       = paramf1.L;
c       = 1; ck = c;    % initial value (c0)
lambda  = 1/L;          % stepsize
N       = 5;           % number of iterations

% store iterates in cells:
x  = cell(N+1,1);  x{1} = x0; 
z  = cell(N+1,1);  z{1} = x0;
y  = cell(N,1); 
g  = cell(N,1); g{1} = gx0;    % gradients of f at the y's
f  = cell(N,1); f{1} = fx0;    % function values of f at the y's
sz = cell(N+1,1); sz{1} = sx0; % subgradients of h at the z's 

for i=1:N
    alphak = (sqrt((ck*lambda)^2+4*ck*lambda)-lambda*ck)/2; 
    ck     = (1-alphak) * ck;
    y{i}   = (1-alphak) * x{i} + alphak * z{i};
    if i>1 % when i = 1, y{i} = x{i} so evaluation is useless
        [g{i},f{i}] = f1.oracle(y{i});
    end
    name    = sprintf('z%d',i);
    z{i+1}  = mirror(g{i}, sz{i}, h, alphak/ck, name); 
    x{i+1} = (1-alphak) * x{i} + alphak * z{i+1};
    
    [sz{i+1}, ~ ] = h.oracle(name); % this is for the next mirror step
end
Dh0 = hxs - hx0 - sx0 * (xs-x0);
P.InitialCondition(c * Dh0 + fx0-fs<=1); % Add an initial condition 

% (4) Set up the performance measure
fN = F.value(x{N+1});         % fN=F(xN)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(fN-fs) 4*L/N^2/c]  % worst-case objective function accuracy
end