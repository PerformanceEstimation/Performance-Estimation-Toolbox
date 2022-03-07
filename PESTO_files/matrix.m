classdef matrix < handle
%MATRIX represents an averaging matrix, used in PESTO for analyzing
%decentralized optimization algorithms.
%Two different matrix description are possible with this class:
%   * Exact description:    It provides the exact expression of the matrix W.
%                           All computations will be exact
%   * Spectral description: It provides a relaxed representation of the
%                           matrix, based on its interval of eigenvalues. 
%                           The represented matrix W is symmetric and doubly stochastic.
%                           The results are a relaxation of PEP.
% Other type of representations (e.g. for non-symmetric matrix) will be
% included in the future.
    
    properties
        type;           % type of matrix description (exact or spectral_relaxed)
        lam;            % spectral matrix description (if it exists)
        time_varying;   % (only for spectral description) can the matrix change over the different consensuses added or not.
        W;              % exact matrix description (if it exists) array of size NxN
        pep;            % instance of PEP problem
        X;              % cell array NxK of PEP variables that are linked to Y using the matrix W (Y=WX)
        Xb;             % cell array 1xK of PEP expression to represents the agent average of variables X, for each consensus K.
        Y;              % cell array NxK of PEP variables that are linked to X using the matrix W (Y=WX)
    end
    
    methods
        function obj = matrix(type,mat,pep,time_varying)
        %MATRIX Construct an instance of this class
        %Input:     - type is either 'exact' or 'spectral_relaxed', depending on the type of description for the matrix
        %           - mat describe the matrix with an NxN array if type='exact', or with a range of eigenvalues (2x1) if type='spectral_relaxed'
        %           - pep is the PEP object of the problem where this matrix will be used.
        %           - time_varying is a boolean that indicates if the different consensus built with this object should use
        %           the same matrix (time_varying = 0), or can use different ones (time_varying = 1)
        %Output:    - the matrix object
            obj.type = type;
            switch type
                case 'spectral_relaxed'
                    assert(length(mat) <= 2, 'for the spectral description of the matrix, mat should be the interval of accepted eigenvalues (2 x 1 array)');
                    assert(all(abs(mat) <= 1), 'for the spectral description of the matrix, the eigenvalues should be between -1 and 1');
                    if length(mat)==1
                        obj.lam = [-mat,mat];
                    else
                        obj.lam = mat;
                    end
                    obj.X = {}; obj.Xb = {}; obj.Y = {};
                    obj.pep = pep;
                    obj.time_varying = time_varying;
                case 'exact'
                    obj.W = mat;
                    obj.pep = pep;
                otherwise
                    assert(0, 'matrix type should be either ''spectral_relaxed'' either ''exact''.');
            end
        end
        
        function y = consensus(obj,x)
        %CONSENSUS creates PEP variables y such that y=Wx,
        % W is a matrix of size NxN (explicitly or spectraly described in this class)
        % INPUT: x is a cell array of size NxK containing PEP variables
        % OUTPUT: y is a cell array of size NxK containing PEP variables
        % This function collects the information about x and creates y
        % but for the spectral description, the spectral constraints imposing that y=Wx
        % should be imposed for all pairs (x,y) at once using the function AddConsensusConstraints.
        % The latter is automatically executed when pep.solve() is called.
            [N, K] = size(x);
            y = cell(N,K);
            switch obj.type
                case 'spectral_relaxed'
                    [sx1,~] = size(obj.X);
                    assert(sx1 == 0 | sx1 == N, 'input x should be a cell array of size NxK, where N is the size of matrix W (NxN) and corresponds to the number of agents');
                    xb = cell(1,K);
                    yb = cell(1,K);
                    for k=1:K
                        for i=1:N   % for each agent
                            % Compute xb
                            if i == 1
                                xb{k} = x{1,k}/N;
                            else
                                xb{k} = xb{k} + x{i,k}/N;
                            end
                            % Creation of Y
                            if i~=N % creation of y variables centered at xb
                                y{i,k} = xb{k} + obj.pep.StartingPoint();
                            elseif N==1
                                y{i,k} = x{i,k};
                            else % setting up last agent such that xb = yb (no new starting points)
                                y{i,k} = N*(xb{k}-yb{k});
                            end
                            % Compute yb (over N-1 agents)
                            if i == 1
                                yb{k} = y{1,k}/N;
                            else
                                yb{k} = yb{k} + y{i,k}/N;
                            end
                        end
                    end
                    obj.X = [obj.X, x]; obj.Xb = [obj.Xb, xb]; obj.Y = [obj.Y, y];
                case 'exact'
                    assert(N == length(obj.W),'input x should be a cell array of size NxK, where N is the size of matrix W (NxN) and corresponds to the number of agents');
                    % compute Y = WX
                    for k=1:K
                        for i=1:N
                            y{i,k} = obj.W(i,1)*x{1,k};
                            for j=2:N
                                y{i,k} = y{i,k} + obj.W(i,j)*x{j,k};
                            end
                        end
                    end
            end
        end
        
        function [Wh,r,status] = estimate(obj, time_varying_mat)
        %ESTIMATE estimates the consensus matrix (or matrices) used for all the
        %consensus steps called on that instance. If time_varying_mat = 1, it
        %return a cell with a different matrix estimates for each consensus
        %step.
        % input:    - time_varying_mat is boolean that indicates if we want
        %           to estimate one unique matrix for all the consensus steps or not.
        %           By default, this argument is set to obj.time_varying.
        % output:   - Wh is the estimate of the consensus matrix such that Y=WX (when time_varying_mat = 0)
        %           when time_varying_mat = 1, Wh is a cell of size 1xK that contains 
        %           a different matrix estimate for each consensus step.
        %           - r is the residual norm of the estimate ||Y-WhX||.
        %           When obj.type='exact', Wh=obj.W and r=0.
        %           When time_varying_mat = 1, r is an array with the different residuals, 
        %           associated with the different matrix estimates, computed on different consensus steps.
        %           - status is a string that gives information about the
        %           status of the estimate.
            if nargin < 2
                time_varying_mat = obj.time_varying;
            end
            switch obj.type 
                case 'spectral_relaxed'
                    [N,K] = size(obj.X);
                    d = length(double(obj.X{1,1}));

                    % Evaluating the X and Y solution of PEP.
                    X_fl = zeros(N,d*K); Y_fl = zeros(N, d*K);
                    for k = 1:K
                        for i = 1:N
                            Y_fl(i,(k-1)*d+1:k*d) = double(obj.Y{i,k});
                            X_fl(i,(k-1)*d+1:k*d) = double(obj.X{i,k});
                        end
                    end
                    if ~time_varying_mat    % constant matrix for each consensus
                        [Wh,r,status] = obj.aux_estimate(X_fl,Y_fl,N);
                    else                    % time-varying matrix
                        Wh = cell(1,K); % different estimates, residuals and status for each consensus
                        r = zeros(1,K);
                        status = cell(1,K);
                        for k = 1:K
                            [Wh{k},r(k),status{k}] = obj.aux_estimate(X_fl(:,(k-1)*d+1:k*d),Y_fl(:,(k-1)*d+1:k*d),N);
                            status{k} = append('time varying: ',status{k});
                        end
                    end
                case 'exact'
                    Wh = obj.W;
                    r = 0;
                    status = 'exact';
            end
        end
        
        function [Wh,r,status] = aux_estimate(obj,X_fl,Y_fl,N)
        %AUX_ESTIMATE is an auxiliary function used in the function ESTIMATE    
            Wh = Y_fl*pinv(X_fl); % pseudo-inverse estimate
            status = "pseudo-inverse estimate";
            % check if Wh is feasible (with tolerance tol)
            tol = 1e-3;
            ev = eig(Wh); ev = ev(~(abs(ev-1)<= tol));
            check_ev = all(obj.lam(2) - ev >= -tol) && all(obj.lam(1) - ev <= tol);
            check_sym = all(abs(Wh - Wh') <= tol,'all');
            check_dstoch = all(abs(ones(1,N)*Wh - ones(1,N)) <= tol);
            if ~check_ev || ~check_sym || ~check_dstoch % Wh not feasible
                 % Compute another estimate based on the following SDP
                 I = eye(N,N); V = sdpvar(N); % symmetric by default
                 objective = norm(Y_fl - V*X_fl,'fro');
                 cons = [];
                 cons = cons + (ones(1,N)*V == ones(1,N)); % stochasticity 
                 cons = cons + (obj.lam(2)*I - (V-1/N*ones(N,N)) >= 0); % bound eigenvalues
                 cons = cons + (obj.lam(1)*I - (V-1/N*ones(N,N)) <= 0); % bound eigenvalues
                 solver_opt = sdpsettings('solver','mosek','verbose',0);
                 optimize(cons,objective,solver_opt);
                 Wh = double(V);
                 status = "SDP estimate";
            end
            r = norm(Wh*X_fl - Y_fl)/sum(size(Y_fl));  % normalized residual norm of the estimate
            if r <= 0.01
                status = status + " - success";
            else
                status = "No worst-matrix";
            end
        end
        
        function AddConsensusConstraints(obj)
        % Note : you shouldn't use this function if you use PESTO normally.
        % Indeed, the function is automatically called inside the function pep.solve()
        % to ensure that it is called before solving the pep problem where matrix W appears.
        % AddConsensusConstraints does nothing when obj.type = 'exact'.
        % When obj.type = 'spectral_relaxed', the consensus matrix is spectrally described and
        % this function adds new constraints to the PEP problem in order to ensure that Y=WX.
        % The new constraints are either scalar inequalities or LMIs depending if the
        % consensus matrix can be time varying or not (obj.time_varying).
            switch obj.type
                case 'spectral_relaxed'
                    [N,K] = size(obj.X);
                    if ~obj.time_varying    % constant matrix for each consensus
                        PSD1 = cell(K,K);
                        PSD2 = cell(K,K);
                        PSD3 = cell(K,K);
                        for k1 = 1:K
                            for k2 = 1:K
                                xx = (obj.X{1,k1}-obj.Xb{k1})*(obj.X{1,k2}-obj.Xb{k2});
                                yy = (obj.Y{1,k1}-obj.Xb{k1})*(obj.Y{1,k2}-obj.Xb{k2});
                                xy = (obj.X{1,k1}-obj.Xb{k1})*(obj.Y{1,k2}-obj.Xb{k2});
                                PSD1{k1,k2} = obj.lam(2)*xx - xy;
                                PSD2{k1,k2} = -obj.lam(1)*xx + xy;
                                PSD3{k1,k2} = -obj.lam(2)*obj.lam(1)*xx + (obj.lam(1)+obj.lam(2))*xy - yy;
                                for i = 2:N
                                    xx = (obj.X{i,k1}-obj.Xb{k1})*(obj.X{i,k2}-obj.Xb{k2});
                                    yy = (obj.Y{i,k1}-obj.Xb{k1})*(obj.Y{i,k2}-obj.Xb{k2});
                                    xy = (obj.X{i,k1}-obj.Xb{k1})*(obj.Y{i,k2}-obj.Xb{k2});
                                    PSD1{k1,k2} = PSD1{k1,k2} + obj.lam(2)*xx - xy;
                                    PSD2{k1,k2} = PSD2{k1,k2} -obj.lam(1)*xx + xy;
                                    PSD3{k1,k2} = PSD3{k1,k2} -obj.lam(2)*obj.lam(1)*xx + (obj.lam(1)+obj.lam(2))*xy - yy;
                                end
                            end
                        end
                        if abs(obj.lam(2)) ~= abs(obj.lam(1))
                            obj.pep.AddLMIConstraint(PSD1); % PSD1 is positive semi-definite
                            obj.pep.AddLMIConstraint(PSD2); % PSD2 is positive semi-definite
                        end
                        obj.pep.AddLMIConstraint(PSD3); % PSD3 is positive semi-definite
                    else                    % time-varying matrix
                        for k=1:K
                            xc2 = (obj.X{1,k} - obj.Xb{k})^2;
                            yc2 = (obj.Y{1,k} - obj.Xb{k})^2;
                            ytx = (obj.Y{1,k}-obj.Xb{k})*(obj.X{1,k}-obj.Xb{k});
                            for i=2:N
                                xc2 = xc2 + (obj.X{i,k}-obj.Xb{k})^2;
                                yc2 = yc2 + (obj.Y{i,k}-obj.Xb{k})^2;
                                ytx = ytx + (obj.Y{i,k}-obj.Xb{k})*(obj.X{i,k}-obj.Xb{k});
                            end
                            if abs(obj.lam(2)) ~= abs(obj.lam(1))
                                obj.pep.AddConstraint(ytx <= obj.lam(2)*xc2);
                                obj.pep.AddConstraint(ytx >= obj.lam(1)*xc2);
                            end
                            obj.pep.AddConstraint(yc2 - (obj.lam(1)+obj.lam(2))*ytx <= -obj.lam(1)*obj.lam(2)*xc2);
                        end
                    end
                otherwise
                  % do nothing      
            end
        end
    end
end