classdef pet < handle
    
    %
    %       Performance Estimation Toolbox (PET) version 20170307
    %
    %       Authors: A. Taylor, J. Hendrickx, F. Glineur
    %       User Manual available at http://perso.uclouvain.be/adrien.taylor
    %       Direct help: help pet
    %       Examples available in the directory pet/examples
    %
    %       Demo's available by typing:
    %       >> demo1 % demo of PET for the gradient descent
    %       >> demo2 % demo of PET for subgradient methods
    %       >> demo3 % demo of PET for FISTA
    %
    
    properties (GetAccess=private)
        expr_list_perf;
        list_size_perf;
        expr_list_init;
        list_size_init;
        expr_list_others;
        list_size_others;
        list_func;
        list_size_func;
        t_reset;
    end
    methods
        function obj=pet()
            obj.expr_list_perf=cell(0,1);
            obj.list_size_perf=0;
            obj.expr_list_init=cell(0,1);
            obj.list_size_init=0;
            obj.expr_list_others=cell(0,1);
            obj.list_size_others=0;
            obj.list_func=cell(0,1);
            obj.list_size_func=0;
            obj.t_reset=now;
            Point.Reset(obj.t_reset);
            pet.CountActive(1);
        end
        function delete(obj)
            pet.CountActive(-1);
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
        function out=AddObjective(obj,InterpEval,param)
            assert(isa(InterpEval,'function_handle') | isa(InterpEval,'char'),'Invalid component added to the objective function');
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
            if isa(InterpEval,'function_handle')
                out=functionClass(InterpEval);
                obj.list_size_func=obj.list_size_func+1;
                obj.list_func{obj.list_size_func,1}=out;
            else
                switch InterpEval
                    case 'Convex'
                        out=functionClass(@(pt1,pt2)Interpolation_Convex(pt1,pt2));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'SmoothStronglyConvex'
                        out=functionClass(@(pt1,pt2)Interpolation_SmoothStronglyConvex(pt1,pt2,param.mu,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'ConvexIndicator'
                        out=functionClass(@(pt1,pt2)Interpolation_ConvexIndicator(pt1,pt2,param.D,param.R));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'ConvexSupport'
                        out=functionClass(@(pt1,pt2)Interpolation_ConvexSupport(pt1,pt2,param.D,param.R));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Smooth'
                        out=functionClass(@(pt1,pt2)Interpolation_Smooth(pt1,pt2,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'StronglyConvexBoundedDomain'
                        out=functionClass(@(pt1,pt2)Interpolation_StronglyConvexBoundedDomain(pt1,pt2,param.D,param.R,param.mu));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'SmoothConvexBoundedGradient'
                        out=functionClass(@(pt1,pt2)Interpolation_SmoothConvexBoundedGradient(pt1,pt2,param.D,param.R,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    otherwise
                        assert(0,'Invalid component added to the objective function');
                end
            end
        end
        function cons=collect(obj,tau,verbose_pet)
            cons=[];
            if obj.list_size_perf>0
                if verbose_pet, fprintf(' PET: Setting up the problem: performance measure'), end;
                for i=1:obj.list_size_perf
                    lexpr=obj.expr_list_perf{i,1};
                    cons=cons+(tau<=lexpr.Eval());
                end
                perf_size=length(cons);
                if verbose_pet, fprintf(' (done, %d constraint(s) added) \n',perf_size), end;
            end
            count=length(cons);
            if obj.list_size_init>0
                if verbose_pet, fprintf(' PET: Setting up the problem: initial conditions'), end;
                for i=1:obj.list_size_init
                    lexpr=obj.expr_list_init{i,1}.Eval();
                    cons=cons+lexpr;
                end
                init_size=length(cons)-count;
                if verbose_pet, fprintf(' (done, %d constraint(s) added) \n',init_size), end;
            end
            count=length(cons);
            if verbose_pet, fprintf(' PET: Setting up the problem: other constraints'), end;
            if obj.list_size_others>0
                for i=1:obj.list_size_others
                    cons=cons+obj.expr_list_others{i,1}.Eval();
                end
            end
            
            for i=1:obj.list_size_func
                cons=cons+obj.list_func{i,1}.collect();
            end
            others_size=length(cons)-count;
            if verbose_pet, fprintf(' (done, %d constraint(s) added) \n',others_size), end;
        end
        function out=solve(obj,verbose_pet,solver_opt)
            
            assert(obj.t_reset==Point.Reset(),'Another PET environment has been initialized, re-generate this PEP before solving it')
            if nargin < 2
                verbose_pet=1;
                verbose_solver=0;
                solver_opt = sdpsettings('verbose',verbose_solver);%,'solver','mosek','mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-10);
            elseif nargin < 3
                verbose_solver=0;
                solver_opt = sdpsettings('verbose',verbose_solver);%,'solver','mosek','mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-10);
            end
            
            dim1=Point.GetSize('Point');
            if verbose_pet, fprintf(' PET: Setting up the problem, size of the main PSD matrix: %d x %d\n',dim1,dim1), end
            G=sdpvar(dim1);
            dim2=Point.GetSize('Function value');
            F=sdpvar(dim2,1);
            tau=sdpvar(1,1);
            obj_func=tau;
            cons=(G>=0);
            Evaluable.SetGetFunc(F);Evaluable.SetGetGram(G);
            msg = sprintf(' PET: Setting up the problem: interpolation constraints (component %d out of %d done)\n', 0,obj.list_size_func);
            if verbose_pet, fprintf(msg), end;
            prevlength = numel(msg);
            for i=1:obj.list_size_func
                cons=cons+obj.list_func{i,1}.GetInterp();
                if verbose_pet
                    msg = sprintf(' PET: Setting up the problem: interpolation constraints (component %d out of %d done)\n', i,obj.list_size_func);
                    fprintf(repmat('\b', 1, prevlength));
                    fprintf(msg);
                    prevlength = numel(msg);
                end
            end
            interp_cons=length(cons)-1;
            if verbose_pet, fprintf('      Total interpolation constraints:  %d \n',interp_cons), end;
            
            cons=cons+obj.collect(tau,verbose_pet);
            
            
            if verbose_pet, fprintf(' PET: Calling SDP solver\n'), end;
            out.solverDetails=optimize(cons,-obj_func,solver_opt);
            out.WCperformance=double(obj_func);
            
            if verbose_pet, fprintf(' PET: Solver output: %7.5e, solution status: %s\n',out.WCperformance,out.solverDetails.info), end;
            if verbose_pet, fprintf(' PET: Post-processing\n'), end;
            
            % Approximating P=[x0 ... gN] from G using Cholesky
            % Decomposition
            [V,D]=eig(double(G));%
            tol=1e-9;%We throw away eigenvalues smaller that 1e-9
            eigenV=diag(D); eigenV(eigenV < tol)=0;
            new_D=diag(eigenV); [~,P]=qr(sqrt(new_D)*V.');
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