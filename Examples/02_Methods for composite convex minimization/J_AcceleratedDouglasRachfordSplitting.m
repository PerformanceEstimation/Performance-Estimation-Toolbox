clear all; clc;
% In this example, we use a fast Douglas-Rachford splitting 
% method for solving the composite convex minimization problem
%   min_x { F(x) = f_1(x) + f_2(x) }
%   (for notational convenience we denote xs=argmin_x F(x);
% where f_2 is L-smooth and f_1 is closed, proper and convex.
%
% The method below is due to
% [1] Panagiotis Patrinos, Lorenzo Stella, and Alberto Bemporad. 
%     "Douglas-Rachford splitting: Complexity estimates and accelerated
%     variants." In 53rd IEEE Conference on Decision and Control (2014)
% where the theory is available for quadratics.
%
% Our notations for the algorithms are as follows:
%
%       x_k     = prox_{\alpha f2}(u_k)
%       y_k     = prox_{\alpha f1}(2*x_k-u_k)
%       w_{k+1} = u_k + \theta (y_k - x_k)
%       if k > 1
%           u{k+1} = w{k+1} + (k-2)/(k+1) * (w{k+1} - w{k});
%       else
%           u{k+1} = w{k+1};
%       end
%
% In Douglas-Rachford schemes, w_{k} converge to a point ws such that
%       xs = prox_{\alpha}(ws) is an optimal point.
% Hence we show how to compute the worst-case behavior of F(y{N})-F(xs)
% given that ||w0 - ws || <= 1.


% (0) Initialize an empty PEP
P = pep();


% (1) Set up the objective function
mf = .1; paramf2.mu = mf;  % Strong convexity parameter
Lf = 1;  paramf2.L  = Lf;  % Smoothness parameter
f1 = P.DeclareFunction('Convex');
f2 = P.DeclareFunction('SmoothStronglyConvex',paramf2);
F  = f1 + f2; % F is the objective function

% (2) Set up the starting point and initial condition
w0      = P.StartingPoint();     % x0 is some starting point
[xs,Fs] = F.OptimalPoint('opt'); % xs is an optimal point, and fs=F(xs)
[g1s,~] = f1.oracle('opt');      % this is f_1'(xs)
[g2s,~] = f2.oracle('opt');      % this is f_2'(xs)

% (3) Algorithm
N = 10;            % number of iterations

alpha = .9;                        % <1/Lf (blindly follow theory of [1])
theta = (1-alpha*Lf)/(1+alpha*Lf); % (blindly follow theory of [1])
ws    = xs + alpha * g2s;          % xs = prox_{\alpha f2} (ws)

P.InitialCondition( (w0-ws)^2 <= 1); % Initial condition ||w0-ws||^2 <= 1

w     = cell(N+1,1);
u     = cell(N+1,1);
x     = cell(N,1);
y     = cell(N,1);
fyval = cell(N,1);
w{1}  = w0;
u{1}  = w0;
for i = 1:N
    x{i}                = proximal_step(u{i},f2,alpha);
    [y{i},~,fyval{i}]   = proximal_step(2*x{i}-u{i},f1,alpha);
    w{i+1}              = u{i} + theta * (y{i}-x{i});
    if i > 1
        u{i+1} = w{i+1} + (i-2)/(i+1) * (w{i+1} - w{i});
    else
        u{i+1} = w{i+1};
    end
end


% (4) Set up the performance measure
F_final = f2.value(y{N})+fyval{N};
P.PerformanceMetric(F_final-Fs);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(F_final-Fs)   

% Theory from [1] on quadratics predicts upper bound 2/(alpha*theta*(N+1)^2)
% 
% For convenience, here are a few outputs from the toolbox
% N_list       = [1:10 15 20 25 30];
% pesto_output = [0.2027 0.1929 0.1839 0.1737 0.1627  0.1514...
%     0.1400 0.1289 0.1182  0.1081 0.0675 0.0413 0.0250 0.0147];
% close all; loglog(N_list,pesto_output,'-r'); hold on;
% loglog(N_list,2./(alpha*theta*(N_list+1).^2),'--b')








