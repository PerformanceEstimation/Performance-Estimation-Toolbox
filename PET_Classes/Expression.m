classdef Expression < Evaluable
    properties (GetAccess=private)
        expression1; % containing its dependence??
        coef1;
        expression2;
        coef2;
        pcoef;
        expr_saved;
        when_saved;
    end
    methods
        function obj=Expression(expr1,coef1,expr2,coef2,pcoef)
            assert(isa(expr1,'Evaluable') && isa(expr2,'Evaluable'),'Must be expression objects (PET class: Expression)');
            assert(strcmp(expr1.getType(),expr2.getType()),'Wrong type combination - incompatible for sum expression (PET class: Expression)');
            
            if nargin < 5
                obj.pcoef=0;
            else
                assert(pcoef==0 | strcmp(expr1.getType(),'Function value'),'Wrong type combination (PET class: PrExpression)');
                obj.pcoef=pcoef;
            end
            obj.type=expr1.getType();
            obj.expression1=expr1;
            obj.coef1=coef1;
            obj.expression2=expr2;
            obj.coef2=coef2;
            obj.expr_saved=[];
            obj.when_saved=0;
        end
        function vec=Eval(obj)
            if isempty(obj.expr_saved) || (obj.when_saved<=Evaluable.update(0))
                if strcmp(obj.type,'Function value')
                    vec=obj.expression1.Eval()*obj.coef1+obj.expression2.Eval()*obj.coef2+obj.pcoef;
                else
                    vec=obj.expression1.Eval()*obj.coef1+obj.expression2.Eval()*obj.coef2;
                end
                obj.expr_saved=vec;
                obj.when_saved=now;
            else
                vec=obj.expr_saved;
            end
        end
    end
end