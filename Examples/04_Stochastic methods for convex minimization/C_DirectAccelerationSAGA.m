function C_DirectAccelerationSAGA
% In this example, we use "SAGA with Sampled Negative Momentum (SSNM)"
% for solving the finite sum minimization problem
%   min_x {F(x)= 1/n (f1(x)+ ... + fn(x)) + h(x) }
% for notational convenience we denote xs=argmin_x F(x).
% Functions f1,...,fn are assumed L-smooth and convex, and h is closed
% proper and mu-strongly convex with a proximal operator available.
%
% This code compute the exact rate for the Lyapunov function from:
% [1] Kaiwen Zhou, Qinghua Ding, Fanhua Shang, James Cheng, Danli Li,
%     Zhi-Quan Luo. "Direct Acceleration of SAGA using Sampled Negative
%     Momentum." (2019)
% (Theorem 1) 

% (0) Initialize an empty PEP
P = pep();

% Problem setup:
L = 1;      % Smoothness
m = .01;    % Strong convexity
n = 5;      % Number of functions in the finite sum

kappa = L/m;

if n/kappa <= 3/4
    eta = 1/sqrt(3*m*L*n);
else
    eta = 2/m/n;
end
tau = n*eta*m/(1+eta*m);


% (1) Set up the objective function
paramf.L   = L;
paramf.mu  = 0;
paramh.L   = Inf;
paramh.mu  = m;

h = P.DeclareFunction('SmoothStronglyConvex',paramh);
F = h;              % F is the objective function
for i = 1:n
    f_i{i} = P.DeclareFunction('SmoothStronglyConvex',paramf);
    F_i{i} = h + f_i{i};
    F     = F + f_i{i}/n;
end
[xs, ~] = F.OptimalPoint('opt');




% (2) Initial values

x{1} = P.StartingPoint();
phi = cell(n,1);
for i = 1:n
    phi{i} = P.StartingPoint();
end

% (3) Initial value of the Lyapunov function

c1 = 1/2/eta/n;
c2 = 1/m/eta/n;
T0 = c1 * (xs-x{1})^2;

for i = 1:n
    name = sprintf('phi%d',i);
    [gi,   fi] = F_i{i}.oracle(phi{i},name);
    [gis, fis] = F_i{i}.oracle('opt');
    T0  = T0 + c2/n * (fi - fis - gis*(phi{i}-xs));
end
P.InitialCondition( T0 <= 1);


% (4) Explicit computation of the expected value of the Lyapunov function 
% after one iteration (so: expectation over n^2 possible scenarios)

T1avg = 0;
for i = 1:n
    name          = sprintf('phi%d',i);
    y{i}          = tau * x{1} + (1-tau) * phi{i};
    grad_estim{i} = f_i{i}.gradient(y{i}) - f_i{i}.gradient(name);
    for j = 1:n
        name = sprintf('phi%d',j);
        grad_estim{i} = grad_estim{i} + 1/n * f_i{j}.gradient(name);
    end
    name          = sprintf('x%d',1+i);
    w{i}          = x{1} - eta * grad_estim{i};
    x{i+1}        = proximal_step( w{i}, h,  eta, name);
    T1            = c1 * (xs-x{1+i})^2;
    for j = 1:n
        for k = 1:n
            new_phi = tau * x{i+1} + (1-tau) * phi{k};
            [gjs, fjs] = F_i{j}.oracle('opt');
            if k ~= j
                name = sprintf('phi%d',j);
                [gj,   fj] = F_i{j}.oracle(name);
                T1  = T1 + 1/n^2 * (fj - fjs - gjs*(phi{j}-xs));
            else
                [gj,   fj] = F_i{k}.oracle(new_phi);
                T1  = T1 + 1/n^2 * (fj - fjs - gjs*(new_phi-xs));
            end
        end
    end
    T1avg = T1avg + T1/n;
end
P.PerformanceMetric(T1avg);


% (5) Solve the PEP
P.solve()

% PESTO output vs Theorem 1 of [1]
[double(T1avg) 1/(1+eta*m)]
end