classdef PrExpression < Evaluable
    properties (GetAccess=private)
        expression1;
        expression2;
        coef;
        pcoef;
        expr_saved;
        when_saved;
    end
    methods
        function obj=PrExpression(expr1,expr2,coef,pcoef)
            assert(isa(expr1,'Evaluable') && isa(expr2,'Evaluable') && isa(coef,'double'),'Must be expression objects (PET class: PrExpression)');
            assert(strcmp(expr1.getType(),expr2.getType()) && strcmp('Point',expr1.getType()),'Wrong type combination - incompatible for product expression (PET class: PrExpression)');
            if nargin < 4
                obj.pcoef=0;
            else
                assert(isa(pcoef,'double'),'Wrong type combination (PET class: PrExpression)');
                obj.pcoef=pcoef;
            end
            obj.type='Function value';%scalar!
            obj.expression1=expr1;
            obj.expression2=expr2;
            obj.coef=coef;
            obj.when_saved=0;
            obj.expr_saved=[];
        end
        function vec=Eval(obj)
            if isempty(obj.expr_saved) || (obj.when_saved<=Evaluable.update(0))
                G=Evaluable.SetGetGram();
                vec=obj.coef*(obj.expression1.Eval().'*G*obj.expression2.Eval())+obj.pcoef;
                obj.expr_saved=vec;
                obj.when_saved=now;
            else
                vec=obj.expr_saved;
            end
        end
    end
end