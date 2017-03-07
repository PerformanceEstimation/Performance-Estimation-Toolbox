classdef pet < handle
    
    %
    %       Performance Estimation Toolbox (PET) version 20170307
    %
    %       Authors: A. Taylor, J. Hendrickx, F. Glineur
    %       User Manual available at http://perso.uclouvain.be/adrien.taylor
    %       Direct help: help pet
    %       Examples available in the directory pet/examples
    %
    %
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
        end
        function disp(obj)
            fprintf('Instance of Performance Estimation Problem (PEP) \n');
            fprintf('Current status: %d component functions, %d initial conditions, %d performance constraints.\n',obj.list_size_func,obj.list_size_init,obj.list_size_perf)
        end
        function pt=GenStartingPoint(obj)
            pt=Point('Point');
        end
        function obj=AddConstraint(obj,expr)
            assert(isa(expr,'Evaluable'),'Invalid initial condition');
            assert(strcmp(expr.getType(),'Point') || strcmp(expr.getType(),'Function value'),'Invalid initial condition');
            obj.list_size_others=obj.list_size_others+1;
            obj.expr_list_others{obj.list_size_others,1}=expr;
        end
        function obj=AddPerformanceConstraint(obj,expr)
            assert(isa(expr,'Evaluable'),'Invalid initial condition');
            assert(strcmp(expr.getType(),'Point') || strcmp(expr.getType(),'Function value'),'Invalid initial condition');
            obj.list_size_perf=obj.list_size_perf+1;
            obj.expr_list_perf{obj.list_size_perf,1}=expr;
        end
        function AddInitialCondition(obj,expr,value)
            assert(isa(expr,'Evaluable') && isa(value,'double'),'Invalid initial condition');
            assert(strcmp(expr.getType(),'Point') || strcmp(expr.getType(),'Function value'),'Invalid initial condition');
            obj.list_size_init=obj.list_size_init+1;
            obj.expr_list_init{obj.list_size_init,1}=expr;
            obj.expr_list_init{obj.list_size_init,2}=value;
        end
        function out=AddComponentObjective(obj,InterpEval,param)
            assert(isa(InterpEval,'function_handle') | isa(InterpEval,'char'),'Invalid component added to the objective function');
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
                        out=functionClass(@(pt1,pt2)Interpolation_ConvexIndicator(pt1,pt2,param.M1,param.M2));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'ConvexSupport'
                        out=functionClass(@(pt1,pt2)Interpolation_ConvexSupport(pt1,pt2,param.M1,param.M2));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'Smooth'
                        out=functionClass(@(pt1,pt2)Interpolation_Smooth(pt1,pt2,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'StronglyConvexBoundedDomain'
                        out=functionClass(@(pt1,pt2)Interpolation_StronglyConvexBoundedDomain(pt1,pt2,param.M1,param.M2,param.mu));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    case 'SmoothConvexBoundedGradient'
                        out=functionClass(@(pt1,pt2)Interpolation_SmoothConvexBoundedGradient(pt1,pt2,param.M1,param.M2,param.L));
                        obj.list_size_func=obj.list_size_func+1;
                        obj.list_func{obj.list_size_func,1}=out;
                    otherwise
                        assert(0,'Invalid component added to the objective function');
                end
            end
        end
        function cons=collect(obj,tau)
            cons=[];
            if obj.list_size_perf>0
                fprintf(' PET: Setting up the problem: performance measure');
                for i=1:obj.list_size_perf
                    lexpr=ExpressionWrapper(obj.expr_list_perf{i,1}.Eval(),obj.expr_list_perf{i,1}.getType());
                    if ~strcmp(lexpr.getType(),'Scalar')
                        lexpr2=lexpr^2;
                        cons=cons+(tau<=lexpr2.Eval());
                    else
                        cons=cons+(tau<=lexpr.Eval());
                    end
                end
                perf_size=length(cons);
                fprintf(' (done, %d constraint(s) added) \n',perf_size);
            end
            count=length(cons);
            if obj.list_size_init>0
                fprintf(' PET: Setting up the problem: initial conditions');
                for i=1:obj.list_size_init
                    lexpr=ExpressionWrapper(obj.expr_list_init{i,1}.Eval(),obj.expr_list_init{i,1}.getType());
                    if ~strcmp(lexpr.getType(),'Scalar')
                        cons=cons+(lexpr^2-obj.expr_list_init{i,2}<=0);
                    else
                        cons=cons+(lexpr-obj.expr_list_init{i,2}<=0);
                    end
                end
                init_size=length(cons)-count;
                fprintf(' (done, %d constraint(s) added) \n',init_size);
            end
            count=length(cons);
            fprintf(' PET: Setting up the problem: other constraints');
            if obj.list_size_others>0
                for i=1:obj.list_size_others
                    lexpr=ExpressionWrapper(obj.expr_list_others{i,1}.Eval(),obj.expr_list_others{i,1}.getType());
                    cons=cons+(lexpr<=0);
                end
            end
            
            for i=1:obj.list_size_func
                cons=cons+obj.list_func{i,1}.collect();
            end
            others_size=length(cons)-count;
            fprintf(' (done, %d constraint(s) added) \n',others_size);
        end
        function out=solve(obj,solver_opt)
            
            if nargin < 2
                verbose=0;
                solver_opt = sdpsettings('verbose',verbose);%,'solver','mosek','mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-10);
            end
            
            dim1=Point.GetSize('Point');
            fprintf(' PET: Setting up the problem, size of the main PSD matrix: %d x %d\n',dim1,dim1);
            G=sdpvar(dim1);
            dim2=Point.GetSize('Function value');
            F=sdpvar(dim2,1);
            tau=sdpvar(1,1);
            obj_func=tau;
            cons=(G>=0);
            ExpressionWrapper.SetGetFunc(F);ExpressionWrapper.SetGetGram(G);
            msg = sprintf(' PET: Setting up the problem: interpolation constraints (component %d out of %d done)\n', 0,obj.list_size_func);
            fprintf(msg);
            prevlength = numel(msg);
            for i=1:obj.list_size_func
                cons=cons+obj.list_func{i,1}.GetInterp();
                msg = sprintf(' PET: Setting up the problem: interpolation constraints (component %d out of %d done)\n', i,obj.list_size_func);
                fprintf(repmat('\b', 1, prevlength));
                fprintf(msg);
                prevlength = numel(msg);
            end
            interp_cons=length(cons)-1;
            fprintf('      Total interpolation constraints:  %d \n',interp_cons);
            
            cons=cons+obj.collect(tau);
            
            
            fprintf(' PET: Calling SDP solver\n');
            out.solverDetails=optimize(cons,-obj_func,solver_opt);
            out.WCperformance=double(obj_func);
            fprintf(' PET: Solver output: %7.5e, solution status: %s\n',out.WCperformance,out.solverDetails.info);
        end
    end
end