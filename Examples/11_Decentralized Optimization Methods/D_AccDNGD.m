function D_AccDNGD
% In this example, we consider K iterations of the Accelerated Distributed 
% Nesterov Gadrient Descent (AccDNGD) [1] with N agents, that each holds a
% local L-smooth and convex function Fi, for solving the following
% decentralized problem: 
%   min_x F(x);     where F(x) is the average of local functions Fi.
% For notational convenience we denote xs=argmin_x F(x), and Fs = F(xs).
% Agents start all with x0 = v0 = y0 = 0, and s0 = avg_i grad F(0).
% Additionnaly, we impose that: avg_i ||xi0 - xs||^2 <= 2
%
% This example shows how to obtain the worst-case performance of AccDNGD
% with PESTO in that case. The performance criterion used is the
% functionnal error at the average of the last iterate: F(avg_i(xiK)) - Fs. 
% Communications between the agents can be formulated in two different ways
% in the PEP problem, leading to different types of worst-case solution:
%   (a) Using a fixed communication network, which leads to exact worst-case results 
%   that are specific to the choosen matrix.
%   (b) Using an entire spectral class of communication networks: the ones represented by 
%   a symmetric (generealized) doubly stochastic matrix with a given range of eigenvalues
%   This leads to a relaxation of PEP, providing worst-case valid for any 
%   matrix of the spectral class.
% Both formulations can be tested here but (b) is used by default.

% For details, see
%   [1] G. Qu and N. Li, “Accelerated distributed nesterov gradient descent”,
%   IEEE Transactions on Automatic Control, 2020.
%   [2] Colla, Sebastien, and Julien M. Hendrickx. "Automatic Performance Estimation
%   for Decentralized Optimization" (2022).

verbose = 1;                % Print the problem set up and the results

%%% Set up general problem parameters %%%
% The system
N = 2;                      % Number of agents

% (a) Exact formulation (fixed network W)
% type = 'exact';
% mat = [0.25,0.75;0.75,0.25]; 
% lam = max(abs(eig(mat-1/N*ones(N,N))));

% (b) Spectral formulation
type = 'spectral_relaxed';  % type of representation for the averaging matrix
lam = 0.5;
mat = [-lam,lam];           % Range of eigenvalues for the symmetric(generalized) doubly stochastic averaging matrix W

% The algorithm
K = 10;                     % Number of iterations of DGD
beta = 0.5;                 % rate of decrease for the step-size
eta = 0.5; k0 = 1; L = 1;        
eta_t = eta./((0:K)+k0).^beta; % diminishing step-size
alpha = sqrt(eta_t(1)*L)*ones(K+1,1); % other parameter of AccDNGD algorithm
equalStart = 1;             % initial iterates are equal for each agent
D = sqrt(2);                % Constants for the initial conditions
time_varying_mat = 0;       % The same averaging matrix is used at each iteration (if 1, no constraints for imposing that the same matrix is used at each iteration)

% (0) Initialize an empty PEP
P = pep();   

% (1) Set up the local and global objective functions
fctClass = 'SmoothStronglyConvex'; % Class of functions to consider for the worst-case
fctParam.L  = L;
fctParam.mu = 0;
returnOpt = 0;
[Fi,Fav,~,~] = P.DeclareMultiFunctions(fctClass,fctParam,N,returnOpt);
[xs,Fs] = Fav.OptimalPoint(); 

% Iterates cells
X  = cell(N,K+1);           % local iterates
V  = cell(N,K+1);
Y  = cell(N,K+1);           
S = cell(N, K);             % S contains the local estimates of the global gradient
F_saved = cell(N,K);
G_saved = cell(N,K);

% (2) Set up the starting points and initial conditions
% each agent start with identical x0 = v0 = y0.
X(:,1) = P.MultiStartingPoints(N,equalStart);
V(:,1) = X(:,1); 
Y(:,1) = X(:,1);
%P.AddConstraint(X{1,1}^2 == 0); This is used in theoretical developments in [1] but it does not impact the convergence results.
% We remove the constraint for numerical stability reasons.

[G_saved(:,1),F_saved(:,1)] = LocalOracles(Fi,X(:,1));
S(:,1) = {1/N*sumcell(G_saved(:,1))};   % s0 = avg_i grad Fi(x0)
P.AddConstraint(1/N*sumcell(foreach(@(xi) (xi-xs)^2,X(:,1))) <= D^2); % avg_i ||xi0 - xs||^2 <= D^2

% (3) Set up the averaging matrix
W = P.DeclareConsensusMatrix(type,mat,time_varying_mat);

% (4) Algorithm (AccDNGD)
% Iterations
for k = 1:K
    c = eta_t(k+1)/eta_t(k)*alpha(k)^2;
    alpha(k+1) = (-c + sqrt(c^2+4*c))/2;
    X(:,k+1) = foreach( @(Wy, S) Wy - eta_t(k)*S, W.consensus(Y(:,k)), S(:,k));
    V(:,k+1) = foreach( @(Wv, S) Wv - eta_t(k)/alpha(k)*S, W.consensus(V(:,k)), S(:,k));
    Y(:,k+1) = foreach( @(X, V) (1-alpha(k+1))*X + alpha(k+1)*V, X(:,k+1), V(:,k+1));
    if k<K
        [G_saved(:,k+1),F_saved(:,k+1)] = LocalOracles(Fi,Y(:,k+1));
        S(:,k+1) = foreach(@(Ws,G2,G1) Ws + G2-G1, W.consensus(S(:,k)), G_saved(:,k+1), G_saved(:,k));
    end
end

% (5) Set up the performance measure
Xav = 1/N*sumcell(X(:,K+1));
[~,FavK] = Fav.oracle(Xav);
P.PerformanceMetric(FavK - Fs); % F(avg_i(xiK)) - Fs

% Activate the trace heuristic for trying to reduce the solution dimension
%P.TraceHeuristic(1);

% (6) Solve the PEP
if verbose
    switch type
        case 'spectral_relaxed'
            fprintf("Spectral PEP formulation for AccDNGD after %d iterations, with %d agents \n",K,N);
            fprintf("Using the following spectral range for the averaging matrix: [%1.2f, %1.2f] \n",mat)
        case 'exact'
            fprintf("Exact PEP formulation for AccDNGD after %d iterations, with %d agents \n",K,N);
            fprintf("The used averaging matrix is\n")
            disp(mat);
    end
end
out = P.solve(verbose+1);
if verbose, out, end

% (7) Evaluate the output
wc = out.WCperformance;

% (8) Construct an approximation of the worst averaging matrix used
[Wh.W,Wh.r,Wh.status] = W.estimate(0);

if verbose
    fprintf("Performance guarantee from PESTO: %1.4f \n\n",wc);
end

end

