classdef Constraint < handle
    properties (GetAccess=private)
        expr;
        signtype;
        next;
    end
    methods
        function obj=Constraint(expr,signtype)
            assert(isa(expr,'Evaluable'),'Constraints must involve scalar values only (PESTO class: Constraint)');
            assert(strcmp(expr.getType(),'Function value'),'Constraints must involve scalar values only (PESTO class: Constraint)');
            assert(strcmp(signtype,'le') || strcmp(signtype,'ge') || strcmp(signtype,'eq'), 'This is not a valid constraint (PESTO class: Constraint)');
            obj.signtype=signtype;
            obj.expr=expr;
            obj.next=[];
        end
        function cons=Eval(obj)
            exprm=obj.expr.Eval();
            cons=[];
            if ~isa(exprm,'double')
                switch obj.signtype
                    case 'le'
                        cons=(exprm<=0);
                    case 'ge'
                        cons=(exprm>=0);
                    case 'eq'
                        cons=(exprm==0);
                end
            else
                if (exprm~=0 && strcmp(obj.signtype,'eq')) || (exprm>0 && strcmp(obj.signtype,'le')) || (exprm<0 && strcmp(obj.signtype,'ge'))
                    assert(false,'Warning, some constraints should be checked (infeasible with no dependence in the variables)');
                end
            end
            if ~isempty(obj.next)
                cons=cons+obj.next.Eval();
            end
        end
        function cons=plus(a,b)
            assert(isa(a,'Constraint') & isa(b,'Constraint'),'Concatenation must involve constraints only (PESTO class: Constraint)');
            cons=a;
            a.next=b;
        end
        function disp(obj)
            fprintf('PESTO Constraint\n');
        end
    end
end