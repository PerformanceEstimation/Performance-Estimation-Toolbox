classdef pep < handle
    
    %
    %       Performance Estimation Toolbox (PESTO) version 20191130
    %
    %       Authors: A. Taylor, J. Hendrickx, F. Glineur
    %       User Manual available as UserGuide.pdf
    %
    %       Direct help: help pesto, help pep
    %       Examples available in the directory Examples/
    %
    %       Demos are available by typing:
    %       >> demo
    %
    %       Reference paper available in PESTO_CDC2017_FINAL.pdf:
    %       Taylor, Adrien B., Julien M. Hendrickx, and Francois Glineur.
    %       "Performance Estimation Toolbox (PESTO): automated worst-case
    %       analysis of first-order optimization methods." Proceedings of
    %       the 56th IEEE Conference on Decision and Control (CDC 2017).
    %
    %   Content:
    %       - Classes of functions: Convex, ConvexIndicator, ConvexSupport,
    %       Smooth, SmoothConvexBoundedGradient, SmoothStronglyConvex,
    %       StronglyConvexBoundedDomain, ConvexBoundedGradient.
    %           Help is available for all functional classes (example: type
    %           'help Convex')
    %       - Primitive steps available: exactlinesearch_step,
    %       gradient_step, linearoptimization_step, projection_step,
    %       proximal_step.
    %           Help is available for all primitive steps (example: type
    %           'help linearoptimization_step')
    %       - Primitive oracles available: inexactsubgradient,
    %       subgradient.
    %           Help is available for all primitive evaluations (example:
    %           type 'help inexactsubgradient')
    %
    
    properties (GetAccess=private)
        expr_list_perf;
        list_size_perf;
        expr_list_init;
        list_size_init;
        expr_list_others;
        list_size_others;
        expr_list_tab_SDP;
        list_size_tab_SDP;
        list_func;
        list_size_func;
        t_reset;
        pesto_opt_traceheuristic;
    end
    methods
        function obj=pep()
            obj.expr_list_perf=cell(0,1);
            obj.list_size_perf=0;
            obj.expr_list_init=cell(0,1);
            obj.list_size_init=0;
            obj.expr_list_others=cell(0,1);
            obj.list_size_others=0;
            obj.expr_list_tab_SDP=cell(0,1);
            obj.list_size_tab_SDP=0;
            obj.list_func=cell(0,1);
            obj.list_size_func=0;
            obj.t_reset=now;
            obj.pesto_opt_traceheuristic = 0;
            Point.Reset(obj.t_reset);
            pep.CountActive(1);
        end
        function delete(obj)
            pep.CountActive(-1);
        end
        function obj=TraceHeuristic(obj,active)
            assert(active==0 || active==1,'Trace heuristic option should be 0/1');
            obj.pesto_opt_traceheuristic = active;
        end
        function disp(obj)
            fprintf('Instance of Performance Estimation Problem (PEP) \n');
            fprintf('Current status: %d component functions, %d initial conditions, %d performance constraints.\n',obj.list_size_func,obj.list_size_init,obj.list_size_perf)
        end
        function pt=StartingPoint(obj)
            pt=Point('Point');
        end
        function pt=GenStartingPoint(obj)
            fprintf('GenStartingPoint is deprecated, consider using StartingPoint instead\n');
            pt=obj.StartingPoint();
        end
        function obj=AddConstraint(obj,expr)
            assert(isa(expr,'Constraint'),'Invalid constraint');
            obj.list_size_others=obj.list_size_others+1;
            obj.expr_list_others{obj.list_size_others,1}=expr;
        end
        function obj=AddSDPConstraint(obj,expr)
        % AddSDPConstraint allows to add an arbitrary SDP constraint to the PEP problem.
        % Input: expr is cell of PEP expression that represents the matrix
        % that should be constraint to be positive semi-definite.
            assert(iscell(expr),'SDP constraint should be given as a cell of expressions');
            obj.list_size_tab_SDP=obj.list_size_tab_SDP+1;
            obj.expr_list_tab_SDP{obj.list_size_tab_SDP,1} = expr;
        end
        function obj=PerformanceMetric(obj,expr)
            assert(isa(expr,'Evaluable'),'Perfomance measures should be scalar values');
            assert(strcmp(expr.getType(),'Function value'),'Perfomance measures should be scalar values');
            obj.list_size_perf=obj.list_size_perf+1;
            obj.expr_list_perf{obj.list_size_perf,1}=expr;
        end
        function obj=AddPerformanceConstraint(obj,expr)
            fprintf('AddPerformanceConstraint is deprecated, consider using PerformanceMetric instead\n');
            obj.PerformanceMetric(expr);
        end
        function InitialCondition(obj,expr)
            assert(isa(expr,'Constraint'),'Invalid constraint');
            obj.list_size_init=obj.list_size_init+1;
            obj.expr_list_init{obj.list_size_init,1}=expr;
        end
        function AddInitialCondition(obj,expr)
            fprintf('AddInitialCondition is deprecated, consider using InitialCondition instead\n');
            obj.InitialCondition(expr);
        end
        function out=AddComponentObjective(obj,InterpEval,param)
            fprintf('AddComponentObjective is deprecated, consider using AddObjective instead\n');
            out=obj.AddObjective(InterpEval,param);
        end
        function out=DeclareFunction(obj,InterpEval,param)
            if nargin == 3
                out=obj.AddObjective(InterpEval,param);
            else
                out=obj.AddObjective(InterpEval);
            end
        end
        function out=AddObjective(obj,InterpEval,param)
            assert(isa(InterpEval,'function_handle') | isa(InterpEval,'char'),'Invalid component added to the objective function');
            if isa(InterpEval,'function_handle')
                out=functionClass(InterpEval);
                obj.list_size_func=obj.list_size_func+1;
                obj.list_func{obj.list_size_func,1}=out;
            else
                if nargin>=3
                    if ~isfield(param,'L')
                        param.L=Inf;
                    end
                    if ~isfield(param,'mu')
                        param.mu=0;
                    end
                    if ~isfield(param,'D')
                        param.D=Inf;
                    end
                    if ~isfield(param,'R')
                        param.R=Inf;
                    end
                else
                    param.L=Inf;
                    param.mu=0;
                    param.D=Inf;
                    param.R=Inf;
                end
                switch InterpEval
                    case 'Convex'
                        out=functionClass(@(pt1,pt2)Convex(pt1,pt2));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Linear'
                        out=functionClass(@(pt1,pt2)Linear(pt1,pt2));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'SmoothStronglyConvex'
                        out=functionClass(@(pt1,pt2)SmoothStronglyConvex(pt1,pt2,param.mu,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'ConvexIndicator'
                        out=functionClass(@(pt1,pt2)ConvexIndicator(pt1,pt2,param.D,param.R));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'ConvexSupport'
                        out=functionClass(@(pt1,pt2)ConvexSupport(pt1,pt2,param.D,param.R));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Smooth'
                        out=functionClass(@(pt1,pt2)Smooth(pt1,pt2,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'StronglyConvexBoundedDomain'
                        out=functionClass(@(pt1,pt2)StronglyConvexBoundedDomain(pt1,pt2,param.D,param.R,param.mu));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'SmoothConvexBoundedGradient'
                        out=functionClass(@(pt1,pt2)SmoothConvexBoundedGradient(pt1,pt2,param.D,param.R,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'ConvexBoundedGradient'
                        out=functionClass(@(pt1,pt2)ConvexBoundedGradient(pt1,pt2,param.D,param.R));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Monotone'
                        out=functionClass(@(pt1,pt2)Monotone(pt1,pt2));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'StronglyMonotone'
                        out=functionClass(@(pt1,pt2)StronglyMonotone(pt1,pt2,param.mu));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Lipschitz'
                        out=functionClass(@(pt1,pt2)Lipschitz(pt1,pt2,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Cocoercive'
                        out=functionClass(@(pt1,pt2)Cocoercive(pt1,pt2,param.beta));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'CocoerciveStronglyMonotone'
                        out=functionClass(@(pt1,pt2)CocoerciveStronglyMonotone(pt1,pt2,param.beta,param.mu));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'LipschitzStronglyMonotone'
                        out=functionClass(@(pt1,pt2)LipschitzStronglyMonotone(pt1,pt2,param.L,param.mu));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    otherwise
                        assert(0,'Invalid component added to the objective function');
                end
            end
        end
        function [cons,names]=collect(obj,tau,verbose_pet)
            cons=[]; names = {};
            if obj.list_size_perf>0
                if verbose_pet>1, fprintf(' PESTO: Setting up the problem: performance measure'), end;
                for i=1:obj.list_size_perf
                    lexpr=obj.expr_list_perf{i,1};
                    cons=cons+(tau<=lexpr.Eval());
                    names{end+1} = sprintf('Perf%d',i);
                end
                perf_size=length(cons);
                if verbose_pet>1, fprintf(' (done, %d constraint(s) added) \n',perf_size), end;
            end
            count=length(cons);
            if obj.list_size_init>0
                if verbose_pet>1, fprintf(' PESTO: Setting up the problem: initial conditions'), end;
                for i=1:obj.list_size_init
                    lexpr=obj.expr_list_init{i,1}.Eval();
                    cons=cons+lexpr;
                    names{end+1} = sprintf('Init%d',i);
                end
                init_size=length(cons)-count;
                if verbose_pet>1, fprintf(' (done, %d constraint(s) added) \n',init_size), end;
            end
            count=length(cons);
            if verbose_pet>1, fprintf(' PESTO: Setting up the problem: other constraints'), end;
            if obj.list_size_others>0
                for i=1:obj.list_size_others
                    cons=cons+obj.expr_list_others{i,1}.Eval();
                    names{end+1} = sprintf('Other%d',i);
                end
                other_size=length(cons)-count;
                if verbose_pet>1, fprintf(' (done, %d constraint(s) added) \n',other_size), end;
            end
            count=length(cons);
            % new SDP constraints
            if obj.list_size_tab_SDP>0
                if verbose_pet>1, fprintf(' PESTO: Setting up the problem: SDP constraints'), end;
                for i=1:obj.list_size_tab_SDP
                    sdpexpr = obj.expr_list_tab_SDP{i,1};
                    len = length(sdpexpr);
                    SDP = sdpvar(len);
                    for i1=1:len
                        for i2=1:len
                            cons = cons + ( SDP(i1,i2) == sdpexpr{i1,i2}.Eval() );
                            names{end+1} = sprintf('SDP%d_build%d',i,(i1-1)*len+i2);
                        end
                    end
                    cons = cons + ( SDP >= 0 );
                    names{end+1} = sprintf('SDP%d',i);
                end
                SDP_size = length(cons)-count;
                if verbose_pet>1, fprintf(' (done, %d constraint(s) added (included %d matrix SDP constraint) \n',SDP_size, obj.list_size_tab_SDP), end;
            end
            for i=1:obj.list_size_func
                cons=cons+obj.list_func{i,1}.collect();
            end
            others_size=length(cons)-count;
            if verbose_pet>1, fprintf(' (done, %d constraint(s) added) \n',others_size), end;
        end
        function out=solve(obj,verbose_pet,solver_opt)
            
            assert(obj.t_reset==Point.Reset(),'Another PEP environment has been initialized, re-generate this PEP before solving it')
            if nargin < 3
                verbose_solver=0;
                solver_opt = sdpsettings('verbose',verbose_solver);
                if nargin < 2
                    verbose_pet = 1;
                end
            end
            
            dim1=Point.GetSize('Point');
            if verbose_pet, fprintf(' PESTO: Setting up the problem, size of the main PSD matrix: %d x %d\n',dim1,dim1), end
            G=sdpvar(dim1);
            dim2=Point.GetSize('Function value');
            F=sdpvar(dim2,1);
            tau=sdpvar(1,1);
            obj_func=tau;
            cons=(G>=0); names = {};
            Evaluable.SetGetFunc(F);Evaluable.SetGetGram(G);
            msg = sprintf(' PESTO: Setting up the problem: interpolation constraints (component %d out of %d done)\n', 0,obj.list_size_func);
            if verbose_pet>1, fprintf(msg), end;
            prevlength = numel(msg);
            for i=1:obj.list_size_func
                [addcons, addnames] = obj.list_func{i,1}.GetInterp();
                cons=cons+addcons; names = [names addnames];
                if verbose_pet>1
                    msg = sprintf(' PESTO: Setting up the problem: interpolation constraints (component %d out of %d done)\n', i,obj.list_size_func);
                    fprintf(repmat('\b', 1, prevlength));
                    fprintf(msg);
                    prevlength = numel(msg);
                end
            end
            interp_cons=length(cons)-1;
            if verbose_pet>1, fprintf('      Total interpolation constraints:  %d \n',interp_cons), end;
            
            [addcons, addnames] = obj.collect(tau,verbose_pet);
            cons=cons+addcons; names = [names addnames];
            
            if verbose_pet>1, fprintf(' PESTO: Calling SDP solver\n'), end;
            
            out.solverDetails=optimize(cons,-obj_func,solver_opt);
            out.WCperformance=double(obj_func);
            
            out.dualvalues = cell(1,length(cons));
            for i=1:length(cons)
                out.dualvalues{i}=dual(cons(i));
            end
            out.dualnames = names;
            
            if verbose_pet, fprintf(' PESTO: Solver output: %7.5e, solution status: %s\n',out.WCperformance,out.solverDetails.info), end; 
            
            if verbose_pet>1, fprintf(' PESTO: Post-processing\n'), end;
            
            if obj.pesto_opt_traceheuristic % to obtain "low-dimensional" wc functions
                if verbose_pet>1, fprintf(' PESTO: Applying resolve with trace heuristic\n'), end;
                cons = cons + (obj_func >= out.WCperformance);
                names = [names 'Fix'];
                optimize(cons,trace(G),solver_opt);
            end
            
            % Approximating P=[x0 ... gN] from G using Cholesky
            % Decomposition
            [V,D]=eig(double(G));%
            tol=1e-5; %Throw away eigenvalues smaller that tol
            eigenV=diag(D); eigenV(eigenV < tol)=0;
            new_D=diag(eigenV); [~,P]=qr(sqrt(new_D)*V.');
            P=P(1:sum(eigenV>0),:);
            Evaluable.Solved(1,double(F),double(P));
        end
    end
    methods (Static)
        function out=CountActive(add)
            persistent nbEval;
            if isempty(nbEval)
                nbEval=0;
            end
            if nbEval+add>=0
                nbEval=nbEval+add;
            end
            out=nbEval;
        end
    end
end

