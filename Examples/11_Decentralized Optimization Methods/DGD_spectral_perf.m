function [wc, Wh, out] = DGD_spectral_perf(K,alpha,N,lam,IC,equalStart,fctClass,fctParam,single_iter_cons,avgAll,verbose)
% DGD_spectral_perf formulates and solves the spectral formulation of PEP for DGD presented in 
%   Colla, Sebastien, and Julien M. Hendrickx. "Automated Worst-Case
%   Performance Analysis of Decentralized Gradient Descent." (2021)
% The function considers the Decentralized Gradient Descent (DGD) after K iterations with N agents
% linked by any symmetric (generalized) doubly-stochastic matrix with a given range of eigenvalues
% The performance measure considered is F(xav)-F(x*), where F is the global function, 
% x* its minimizer and xav an average of different iterates (all or only the last ones)
% INPUTS:
%   K:          Number of iterations of DGD.
%   alpha:      Array 1 x K with the step-size to use for each iteration. 
%               If alpha is a single number, the step-size is considered as constant.
%   N:          Number of agents in DGD.
%   lam:        Array 1 x 2 describing the range of eigenvalues of the communication matrix W used in DGD: 
%               lam(1) <= lam_i(W) <= lam(2) for i=2,...,N (lam_1(W) = 1).
%   IC:         Constant for the initial condition: ||x_i^0 - x*||^2 <= IC.
%   equalStart: Boolean for deciding if all the agents start with the same initial point (1) or not (0).
%   fctClass:   Name of a PEP functional class. For example 'ConvexBoundedGradient'.
%   fctParam:   Struct containing the parameters of the functional class. For example fctParam.R = 1.
%   single_iter_cons: 
%               Boolean for deciding if the constraints of the spectral formulation 
%               should link only one iteration (1) or all the iterations together (0).
%   avgAll:     Boolean for deciding if 'xav' used in the performance measure averages the iterates of all agents,
%               at each iterations (1) or only the iterates of all agents after the last iteration (0).
%   verbose:    Boolean for deciding to print details about the problem and its resolution (1) or not (0).
% OUTPUTS:
%   wc:         Value of the worst-case performance.
%   Wh:         Approximation of the worst communication matrix.
%   out:        output of the resolution of PESTO.

    % PEP problem
    P = pep();
    if length(lam) == 1
        lam = [lam, lam];
    end    
    
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
    [xs,Fs] = Ftot.OptimalPoint(); 		 % xs is an optimal point, and Fs=F(xs)

    % Iterates cells
    X = cell(N, K+1);
    Y = cell(N, K); % Y = WX
    xb = cell(K,1); % average X over the agents at each ietration
    yb = cell(K,1); % average Y over the agents at each ietration
    
    % Set up the starting point and initial condition
    X{1,1} = P.StartingPoint(); 
    P.InitialCondition((X{1,1}-xs)^2 <= IC^2); % Initial condition ||x_i^0-x*||^2<= IC^2

    xav = X{1,1};       % average value of x among all agents and all iterations

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
    for k = 1:K     % for each iteration
        for i=1:N   % for each agent
            % Update
            Y{i,k} = P.StartingPoint();
            X{i,k+1} = Y{i,k} - alpha(k)*F{i}.subgradient(X{i,k});
            
            % Compute xb and yb
            if i == 1
                xb{k} = X{1,k}/N; yb{k} = Y{1,k}/N;
            else
                xb{k} = xb{k} + X{i,k}/N;
                yb{k} = yb{k} + Y{i,k}/N;
            end      
            xav = xav + X{i,k+1};
        end
        
        % Constraints such that Y = WX
        
        P.AddConstraint((xb{k}-yb{k})^2==0);
        
        % Scalar constraints over only 1 iteration
        if single_iter_cons
            xc2 = (X{1,k} - xb{k})^2;
            yc2 = (Y{1,k} - yb{k})^2;
            ytx = (Y{1,k}-yb{k})*(X{1,k}-xb{k});
            for i=2:N
                xc2 = xc2 + (X{i,k}-xb{k})^2;
                yc2 = yc2 + (Y{i,k}-yb{k})^2;
                ytx = ytx + (Y{i,k}-yb{k})*(X{i,k}-xb{k});
            end
            P.AddConstraint(ytx <= lam(2)*xc2); 
            P.AddConstraint(ytx >= lam(1)*xc2);
            P.AddConstraint(yc2 - (lam(1)+lam(2))*ytx <= -lam(1)*lam(2)*xc2);
        end        
    end
    if ~single_iter_cons % LMIs constraints linking all the iterations together
        LMI1 = cell(K,K);
        LMI2 = cell(K,K);
        LMI3 = cell(K,K);
        for k1 = 1:K
            for k2 = 1:K
                xx = (X{1,k1}-xb{k1})*(X{1,k2}-xb{k2});
                yy = (Y{1,k1}-yb{k1})*(Y{1,k2}-yb{k2});
                xy = (X{1,k1}-xb{k1})*(Y{1,k2}-yb{k2});
                LMI1{k1,k2} = lam(2)*xx - xy;
                LMI2{k1,k2} = -lam(1)*xx + xy;
                LMI3{k1,k2} = -lam(2)*lam(1)*xx + (lam(1)+lam(2))*xy - yy;
                for i = 2:N
                    xx = (X{i,k1}-xb{k1})*(X{i,k2}-xb{k2});
                    yy = (Y{i,k1}-yb{k1})*(Y{i,k2}-yb{k2});
                    xy = (X{i,k1}-xb{k1})*(Y{i,k2}-yb{k2});
                    LMI1{k1,k2} = LMI1{k1,k2} + lam(2)*xx - xy;
                    LMI2{k1,k2} = LMI2{k1,k2} -lam(1)*xx + xy;
                    LMI3{k1,k2} = LMI3{k1,k2} -lam(2)*lam(1)*xx + (lam(1)+lam(2))*xy - yy;
                end
            end            
        end   
      P.AddLMIConstraint(LMI1); % LMI1 is positive semi-definite
      P.AddLMIConstraint(LMI2); % LMI2 is positive semi-definite
      P.AddLMIConstraint(LMI3); % LMI3 is positive semi-definite
    end


    % Set up the performance measure
    if avgAll % Average over all agents, all iterations
        xav = xav/((K+1)*N);
        P.PerformanceMetric(Ftot.value(xav)-Fs); % Worst-case evaluated as F(xav)-F(x*)
    else      % Average only for the last iteration
        xav = X{1,K+1}/N;
        for i = 2:N
            xav = xav + X{i,K+1}/N;
        end
        P.PerformanceMetric(Ftot.value(xav)-Fs);
    end
    
    % Activate the trace heuristic for trying to reduce the solution dimension
    %P.TraceHeuristic(1);
    
    % Solve the PEP
    if verbose
        fprintf("Spectral PEP formulation for DGD after %d iterations, with %d agents \n",K,N);
        fprintf("Using the following spectral range for the communication matrix: [%1.2f, %1.2f] \n",lam)
    end
    out = P.solve(verbose);

    % Evaluate the output
    wc = out.WCperformance;
    
    % Construct an approximation of the worst communication matrix that links the solutions X and Y
    d = length(double(X{1,1}));
    % Compute the double value of the iterates
    Xf = zeros(N,d*K); Yf = zeros(N, d*K);
    for k = 1:K
        for i = 1:N
            Xf(i,(k-1)*d+1:k*d) = double(X{i,k});
            Yf(i,(k-1)*d+1:k*d) = double(Y{i,k});
        end
    end
    Wh = Yf*pinv(Xf); % Using pseudo-inverse
end