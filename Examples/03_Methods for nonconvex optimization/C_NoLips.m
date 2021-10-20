function C_NoLips
% In this example, we use a Bregman gradient method for
% solving the constrained nonconvex minimization problem
%   min_x { F(x) = f_1(x) + f_2(x) }
%   for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is h-smooth (possibly nonconvex) and where f_2(x) is a
% closed convex indicator function.
%
% We show how to compute the worst-case value of 
%               min_{1<=i<=N} Dh(x_{i-1},x_i) 
% when xN is obtained by doing N steps of the method and under
% the condition F(x0)-F(xN)<=1. 
%
% [1] Jérôme Bolte, Shoham Sabach, Marc Teboulle, and Yakov Vaisbourd.
%     "First order methods beyond convexity and Lipschitz gradient
%     continuity with applications to quadratic inverse problems."  (2018)
%
%J ́erˆome Bolte∗Shoham Sabach†Marc Teboulle‡Yakov Vaisbourd
% [2] Radu-Alexandru Dragomir, Adrien B. Taylor, Alexandre d’Aspremont, and
%     Jérôme Bolte. "Optimal Complexity and Certification of Bregman
%     First-Order Methods". (2019)
%
%
% DISCLAIMER: This examples requires some experience with PESTO and PEPs
% (see Section 4 in [2]).
%

% (0) Initialize an empty PEP
P = pep();

L   = 1; % d1 = Lh - f1 is convex and  d2 = f1 + Lh is convex
d1  = P.DeclareFunction('Convex');
d2  = P.DeclareFunction('Convex');
f1  = (d2 - d1)/2;
h   = (d1 + d2)/2/L;
f2  = P.DeclareFunction('ConvexIndicator');

F  = f1 + f2;
% (2) Set up the starting point and initial condition
x0        = P.StartingPoint();        % x0 is some starting point

[gh0, h0] = h.oracle(x0,'x0');
[gf0, f0] = f1.oracle(x0,'x0');
[~, F0]   = F.oracle(x0,'x0');

% (3) Algorithm
gamma = 1/L;    % stepsize
N     = 3;     % number of iterations

x   = cell(N+1,1);   x{1} = x0;
gfx = cell(N+1,1); gfx{1} = gf0;
fx  = cell(N+1,1);  fx{1} = f0;
ghx = cell(N+1,1); ghx{1} = gh0;
hx  = cell(N+1,1);  hx{1} = h0;
Dhi = cell(N  ,1);
for i = 1:N
    name    = sprintf('x%d',i);
    x{i+1}	= mirror(gfx{i}, ghx{i}, h+f2, gamma, name);
    
    [gfx{i+1}, fx{i+1}] = f1.oracle(name);
    [ghx{i+1}, hx{i+1}] = h.oracle(name);
    Dhi{i} = hx{i}-hx{i+1}-ghx{i+1}*(x{i}-x{i+1});
    P.PerformanceMetric(Dhi{i});
end
[~,FN] = F.oracle(x{N+1},'xN');
P.InitialCondition(F0-FN <= 1);    % Initial condition Dh(x*,x0)<=1

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
for i = 1:N
    DhValues(i)=double(Dhi{i});
end
[min(DhValues) gamma/N ]
% Result should match gamma/N
end