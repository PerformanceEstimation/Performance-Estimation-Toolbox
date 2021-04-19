function [wc, out] = DGD_exact_perf(K,alpha,N,W,IC,equalStart,fctClass,fctParam,avgAll,verbose)
% DGD_exact_perf formulates and solves the exact PEP problem for the
% Decentralized Gradient Descent (DGD) after K iterations with N agents linked by matrix A.
% The performance measure considered is F(xav)-F(x*), where F is the global function,
% x* its minimizer and xav an average of different iterates (all or only the last ones)
% INPUTS:
%   K:          Number of iterations of DGD.
%   alpha:      Array 1 x K with the step-size to use for each iteration. 
%               If alpha is a single number, the step-size is considered as constant.
%   N:          Number of agents in DGD.
%   W:          Adjacency matrix (N x N) of the communication network between the agent.
%   IC:         Constant for initial condition: ||x_i^0 - x*||^2 <= IC.
%   equalStart: Boolean for deciding if all the agents start with the same initial point (1) or not (0).
%   fctClass:   Name of a PEP functional class. For example 'ConvexBoundedGradient'.
%   fctParam:   Struct containing the parameters of the functional class. For example fctParam.R = 1.
%   avgAll:     Boolean for deciding if 'xav' used in the performance measure averages the iterates of all the agents,
%               at each iterations (1) or only the iterates of all the agents after the last iteration (0).
%   verbose:    Boolean for deciding to print details about the problem and its resolution (1) or not (0).
% OUTPUTS:
%   wc:         Value of the worst-case performance
%   out:        output of the resolution of PESTO.

    % PEP problem
    P = pep();

    % Set up the local objective functions of each agent
    F = cell(1,N);
    for i=1:N
        F{i}=P.DeclareFunction(fctClass,fctParam); % F is the objective function
    end

    % Set up the global function and its optimum
    Ftot = F{1}/N;
    for i=2:N
        Ftot = Ftot + F{i}/N;
    end
    [xs,Fs] = Ftot.OptimalPoint();  % xs is an optimal point, and Fs=F(xs)

    X = cell(N,K+1);     % Iterates cell

    % Set up the starting points and initial conditions
    X{1,1} = P.StartingPoint();
    xav = X{1,1};       % average value of x among all agents and all iterations
    
    P.InitialCondition((X{1,1}-xs)^2 <= IC^2); % Initial condition ||x_i^0-x*||^2<= IC^2; for all i

    for i=2:N
        if equalStart   % All the agents have the same starting point
            X{i,1} = X{1,1}; 
        else            % The agents have different starting points
            X{i,1} = P.StartingPoint();
            P.InitialCondition((X{i,1}-xs)^2 <= IC^2);
        end
        xav = xav + X{i,1};
    end
    
    % Algorithm DGD
    % Step-size 
    if length(alpha) == 1 % use the same step-size for each iteration
        alpha = alpha*ones(1,K);
    end
    % Iterations
    for k=1:K       % for each iteration
       for i=1:N    % for each agent
            X{i,k+1} = - alpha(k)*F{i}.subgradient(X{i,k}); % move along the local subgradient
            for j=1:N % consensus between the agents
               X{i,k+1} = X{i,k+1} + W(i,j)*(X{j,k});
            end
            xav = xav + X{i,k+1};    
       end
    end
    
    % Set up the performance measure
    if avgAll   % Average over all agents, all iterations
        xav = xav/((K+1)*N);
        P.PerformanceMetric(Ftot.value(xav)-Fs); % Worst-case evaluated as F(xav)-F(x*)
    else        % Average only for the last iteration
        xavLast = X{1,K+1}/N;
        for i = 2:N
            xavLast = xavLast + X{i,K+1}/N;
        end
        P.PerformanceMetric(Ftot.value(xavLast)-Fs); 
    end
    
    % Activate the trace heuristic for trying to reduce the solution dimension
    %P.TraceHeuristic(1); 

    % Solve the PEP
    if verbose
        fprintf("Exact PEP formulation for DGD after %d iterations, with %d agents \n",K,N);
        fprintf("The used communication matrix is")
        W
    end
    out = P.solve(verbose);

    % Evaluate the output
    wc = out.WCperformance;
    
end
