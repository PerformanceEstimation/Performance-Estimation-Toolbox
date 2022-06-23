function B_DIGing
% In this example, we consider K iterations of the DIGing algorithm [1]
% with N agents, that each holds a local L-smooth mu-strongly convex
% function Fi, for solving the following decentralized problem:
%   min_x F(x);     where F(x) is the average of local functions Fi.
% For notational convenience we denote xs=argmin_x F(x).
% Each agent i starts with its own iterates xi0, si0 such that 
%   avg_i ||xi0 - xs||^2 <= 1, si0 = grad Fi(xi0) and avg_i ||si0 - avg_i(grad Fi(xi0))||^2 <= 1
%
% This example shows how to obtain the worst-case performance of DIGing with PESTO in that case.
% The performance criterion used is avg_i ||xiK - xs||^2.
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
%   [1] A. Nedic, A. Olshevsky, and W. Shi, “Achieving geometric convergence
%   for distributed optimization over time-varying graphs,” SIAM Journal on
%   Optimization, 2016.
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
lam = 0.9;
mat = [-lam,lam];           % Range of eigenvalues for the symmetric(generalized) doubly stochastic averaging matrix W

% The algorithm
K = 10;                     % Number of iterations of DIGing
alpha = 0.44*(1-lam)^2;     % Step-size used in DIGing (constant) (hand-tuned formula)
equalStart = 0;             % initial iterates are not necessarily equal for each agent
D = 1; E = 1;               % Constants for the initial conditions
time_varying_mat = 0;       % The same averaging matrix is used at each iteration (if 1, no constraints for imposing that the same matrix is used at each iteration)

% (0) Initialize an empty PEP
P = pep();   

% (1) Set up the local and global objective functions
fctClass = 'SmoothStronglyConvex'; % Class of functions to consider for the worst-case
fctParam.L  = 1;
fctParam.mu = 0.1;
returnOpt = 0;
[Fi,Fav,~,~] = P.DeclareMultiFunctions(fctClass,fctParam,N,returnOpt);
[xs,~] = Fav.OptimalPoint(); 

% Iterates cells
X = cell(N, K+1);           % local iterates
S = cell(N, K);             % S contains the local estimates of the global gradient
F_saved = cell(N,K+1);
G_saved = cell(N,K+1);

% (2) Set up the starting points and initial conditions
X(:,1) = P.MultiStartingPoints(N,equalStart);
[G_saved(:,1),F_saved(:,1)] = LocalOracles(Fi,X(:,1));
P.AddConstraint(1/N*sumcell(foreach(@(xi) (xi-xs)^2,X(:,1))) <= D^2); % avg_i ||xi0 - xs||^2 <= D^2
S(:,1) = G_saved(:,1);
P.AddConstraint(1/N*sumcell(foreach(@(si) (si - 1/N*sumcell(G_saved(:,1)))^2,S(:,1))) <= E^2);
                %avg_i ||si0 - avg_i(grad Fi(xi0))||^2 <= E^2
                
% (3) Set up the averaging matrix
W = P.DeclareConsensusMatrix(type,mat,time_varying_mat);

% (4) Algorithm (DIGing)
% Iterations
for k = 1:K
    X(:,k+1) = foreach(@(Wx,S) Wx-alpha*S, W.consensus(X(:,k)), S(:,k));
    [G_saved(:,k+1),F_saved(:,k+1)] = LocalOracles(Fi,X(:,k+1));
    if k<K
        S(:,k+1) = foreach(@(Ws,G2,G1) Ws + G2-G1, W.consensus(S(:,k)), G_saved(:,k+1), G_saved(:,k)); 
    end
end

% (5) Set up the performance measure
metric = 1/N*sumcell(foreach(@(xiK)(xiK-xs)^2,X(:,K+1))); % average squared error: avg_i ||xiK - xs||^2
P.PerformanceMetric(metric); 

% Activate the trace heuristic for trying to reduce the solution dimension
%P.TraceHeuristic(1);

% (6) Solve the PEP
if verbose
    switch type
        case 'spectral_relaxed'
            fprintf("Spectral PEP formulation for DIGing after %d iterations, with %d agents \n",K,N);
            fprintf("Using the following spectral range for the averaging matrix: [%1.2f, %1.2f] \n",mat)
        case 'exact'
            fprintf("Exact PEP formulation for DIGing after %d iterations, with %d agents \n",K,N);
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

% (9) Comparison with theoretical guarantee [1, Theorem 3.14]
% This guarantee always applies to a spectral class of averaging matrices.
wc_theo = max(sqrt(1-alpha*fctParam.mu/1.5),(sqrt(alpha*2*fctParam.L*(1+4*sqrt(N)*sqrt(fctParam.L/fctParam.mu))) + lam))^K;
msg_theo = '';
if wc_theo >= 1
    msg_theo = 'none,';
end
if verbose
    fprintf("Performance guarantee from PESTO: %1.4f \n",wc);
    fprintf("Theoretical guarantee from [1]: %s %1.4f\n\n",msg_theo,wc_theo);
end
end
