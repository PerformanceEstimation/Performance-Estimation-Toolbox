classdef (Abstract) Evaluable < handle
    properties (GetAccess=public)
        type;
        mult_const;
    end
    methods
        function obj=Evaluable()
            Evaluable.update(1);
            Evaluable.CountActive(1);
        end
        function delete(obj)
            Evaluable.CountActive(-1);
        end
        function expr=getType(obj)
            expr=obj.type;
        end
        function disp(obj)
            if (strcmp(obj.getType(),'Point'))
                msg='PEP variable (type: decision variable or gradient)\n';
            else
                msg='PEP variable (type: function value)\n';
            end
            fprintf(msg);
        end
        function c=le(a,b)
            Evaluable.check_for_ineq(a,b)
            c=Constraint(a-b,'le');
        end
        function c=ge(a,b)
            Evaluable.check_for_ineq(a,b)
            c=Constraint(a-b,'ge');
        end
        function c=eq(a,b)
            Evaluable.check_for_ineq(a,b)
            c=Constraint(a-b,'eq');
        end
        function tf=isEqual(a,b)
            tf=min(a.Eval()-b.Eval()==0);
        end
        function obj3=plus(a,b)
            assert((isa(a,'Evaluable') | isa(a,'double')) & (isa(b,'Evaluable') | isa(b,'double')),'Invalid use of SUM - elements are not compatible (PEsTo class: Evaluable)');
            if (isa(a,'Evaluable') && isa(b,'double'))
                assert(strcmp(a.getType(),'Function value'),'Invalid use of SUM - elements are not of the same dimensions (PEsTo class: Evaluable)');
                obj3=Expression(a,1,Point(a.getType(),0),0,b);
            elseif (isa(b,'Evaluable') && isa(a,'double'))
                assert(strcmp(b.getType(),'Function value'),'Invalid use of SUM - elements are not of the same dimensions (PEsTo class: Evaluable)');
                obj3=Expression(b,1,Point(b.getType(),0),0,a);
            else
                assert(strcmp(a.getType(),b.getType()),'Invalid use of SUM - elements are not of the same dimensions (PEsTo class: Evaluable)')
                obj3=Expression(a,1,b,1,0);
            end
        end
        function obj3=uplus(a)
            obj3=a;
        end
        function obj3=minus(a,b)
            assert((isa(a,'Evaluable') | isa(a,'double')) & (isa(b,'Evaluable') | isa(b,'double')),'Invalid use of SUM - elements are not compatible (PEsTo class: Evaluable)');
            if (isa(a,'Evaluable') && isa(b,'double'))
                assert(strcmp(a.getType(),'Function value'),'Invalid use of SUM - elements are not of the same dimensions (PEsTo class: Evaluable)');
                obj3=Expression(a,1,Point(a.getType(),0),0,-b);
            elseif (isa(b,'Evaluable') && isa(a,'double'))
                assert(strcmp(b.getType(),'Function value'),'Invalid use of SUM - elements are not of the same dimensions (PEsTo class: Evaluable)');
                obj3=Expression(b,-1,Point(b.getType(),0),0,a);
            else
                assert(strcmp(a.getType(),b.getType()),'Invalid use of SUM - elements are not of the same dimensions (PEsTo class: Evaluable)')
                obj3=Expression(a,1,b,-1,0);
            end
        end
        function obj3=uminus(a)
            obj3=Expression(a,-1,Point(a.getType(),0),0,0);
        end
        function obj3=mtimes(a,b)
            assert((isa(a,'Evaluable') | isa(a,'double')) & (isa(b,'Evaluable') | isa(b,'double')),'Invalid use of TIMES - elements are not compatible (PEsTo class: Evaluable)');
            if (isa(a,'Evaluable') && isa(b,'double'))
                obj3=Expression(a,b,Point(a.getType(),0),0,0);
            elseif (isa(b,'Evaluable') && isa(a,'double'))
                obj3=Expression(b,a,Point(b.getType(),0),0,0);
            else
                assert(strcmp(a.getType(),b.getType()),'Invalid use of TIMES - elements are not of the same dimensions (PEsTo class: Evaluable)')
                assert(strcmp(a.getType(),'Point'),'Invalid use of TIMES - scalar products are only defined for coordinates/gradients (PEsTo class: Evaluable)')
                obj3=PrExpression(a,b,1,0);
            end
        end
        function obj3=times(a,b)
            obj3=mtimes(a,b);
        end
        function obj3=mpower(a,b)
            assert(b==2,'Only squares are accepted in expressions (PEsTo class: Evaluable)');
            obj3=mtimes(a,a);
        end
        function obj3=power(a,b)
            obj3=mpower(a,b);
        end
        function obj3=mrdivide(a,b)
            assert(isa(b,'double'),'Invalid use of RDIVIDE - elements are not compatible (PEsTo class: Evaluable)');
            obj3=Expression(a,1/b,Point(a.getType(),0),0,0);
        end
        function obj3=rdivide(a,b)
            obj3=mrdivide(a,b);
        end
        function obj3=ctranspose(a)
            obj3=a;
        end
        function obj3=transpose(a)
            obj3=a;
        end
    end
    methods (Abstract)
        double(obj)
        Eval(obj)
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
        function out=update(upd)
            persistent last_update;
            if upd
                last_update=now;
            elseif isempty(last_update)
                last_update=0;
            end
            out=last_update;
        end
        function check_for_ineq(a,b)
            assert((isa(a,'Evaluable') | isa(a,'double')) & (isa(b,'Evaluable') | isa(b,'double')),'Constraints must involve scalar values only (PEsTo class: Constraint)');
            assert(~(isa(a,'double') & isa(b,'double')),'Constraints must involve variables (PEsTo class: Constraint)');
            if isa(a,'double')
                assert(strcmp(b.getType(),'Function value'),'Constraints must involve scalar values only (PEsTo class: Constraint)');
            elseif isa(b,'double')
                assert(strcmp(a.getType(),'Function value'),'Constraints must involve scalar values only (PEsTo class: Constraint)');
            else
                assert(strcmp(a.getType(),'Function value'),'Constraints must involve scalar values only (PEsTo class: Constraint)');
                assert(strcmp(b.getType(),'Function value'),'Constraints must involve scalar values only (PEsTo class: Constraint)');
            end
        end
        function G=SetGetGram(G)
            persistent Gram;
            if nargin == 1
                Gram=G;
            else
                assert(~isempty(Gram),'Warning, the YALMIP variables were not initialized: impossible to evaluate expressions (Eval was probably called outside a PEsTo instance)')
                G=Gram;
            end
        end
        function F=SetGetFunc(F)
            persistent Func;
            if nargin == 1
                Func=F;
            else
                assert(~isempty(Func),'Warning, the YALMIP variables were not initialized: impossible to evaluate expressions (Eval was probably called outside a PEsTo instance)')
                F=Func;
            end
        end
        function out=Solved(in,F,P)
            persistent isSolved;
            if nargin>=1
                isSolved.status=in;
            elseif isempty(isSolved)
                isSolved.status=0;
            end
            if nargin>=3
                isSolved.F=F;
                isSolved.P=P;
            end
            out=isSolved;
        end
    end
    
end

