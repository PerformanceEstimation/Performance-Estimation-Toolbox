function A_SAGA
% In this example, we use SAGA for solving the finite sum minimization
% problem
%   min_x {F(x)= 1/n (f1(x)+ ... + fn(x)) + h(x) }
% for notational convenience we denote xs=argmin_x F(x).
% Functions f1,...,fn are assumed L-smooth and mu-strongly convex, and h
% is closed proper and convex with a proximal operator available.
%
% This code compute the exact rate for the Lyapunov function from the
% original SAGA paper:
% [1] Aaron Defazio, Francis Bach, and Simon Lacoste-Julien.
%     "SAGA: A fast incremental gradient method with support for
%     non-strongly convex composite objectives." (2014)
% (Theorem 1 of [1])

% (0) Initialize an empty PEP
P = pep();

% Problem setup:
L = 1;	% Smoothness
m = .1; % Strong convexity
n = 5;  % Number of functions in the finite sum

gamma = 1/2/(m*n+L);                % stepsize
c     = 1/2/gamma/(1-m*gamma)/n;    % Lyapunov function parameter
kappa = 1/gamma/m;

% (1) Set up the objective function
param.L   = L;
param.mu  = m;

h = P.DeclareFunction('Convex');
F = h;              % F is the objective function
for i = 1:n
    f{i} = P.DeclareFunction('SmoothStronglyConvex',param);
    F    = F + f{i}/n;
end
[xs, ~] = F.OptimalPoint('opt');

% (2) Initial values
phi = cell(n,1);
for i = 1:n
    phi{i} = P.StartingPoint();
end
x{1} = P.StartingPoint();

% (3) Initial value of the Lyapunov function
T0 = c * (xs-x{1})^2;

for i = 1:n
    name = sprintf('phi%d',i);
    [gi,   fi] = f{i}.oracle(phi{i},name);
    [gis, fis] = f{i}.oracle('opt');
    T0  = T0 + 1/n * (fi - fis - gis*(phi{i}-xs));
end

P.InitialCondition( T0 <= 1);

% (4) Explicit computation of the expected value of the Lyapunov function
% after one iteration (so: expectation over n possible scenarios: one for 
% each element fi in the function).

T1avg = 0;
for i = 1:n
    name = sprintf('phi%d',i);
    w{i} = x{1} - gamma * ( f{i}.gradient(x{1},'x0') - f{i}.gradient(name) );
    for j = 1:n
        name = sprintf('phi%d',j);
        w{i} = w{i} - gamma/n * f{j}.gradient(name);
    end
    x{i+1} = proximal_step( w{i}, h,  gamma);
    T1 = c * (xs-x{1+i})^2;
    for j = 1:n
        [gis, fis] = f{j}.oracle('opt');
        if i ~= j
            name = sprintf('phi%d',j);
            [gi,   fi] = f{j}.oracle(name);
            T1  = T1 + 1/n * (fi - fis - gis*(phi{j}-xs));
        else
            name = sprintf('x0');
            [gi,   fi] = f{j}.oracle(name);
            T1  = T1 + 1/n * (fi - fis - gis*(x{1}-xs));
        end
    end
    T1avg = T1avg + T1/n;
end
P.PerformanceMetric(T1avg);

% (5) Solve the PEP
P.solve()

% PESTO output vs Theorem 1 of [1]
[double(T1avg) (1-1/kappa)]
end

