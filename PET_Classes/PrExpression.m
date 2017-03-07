classdef PrExpression < Evaluable
    properties (GetAccess=private)
        expression1;
        expression2;
        coef;
        pcoef;
        type;
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
                G=ExpressionWrapper.SetGetGram();
                vec=obj.coef*(obj.expression1.Eval().'*G*obj.expression2.Eval())+obj.pcoef;
                obj.expr_saved=vec;
                obj.when_saved=now;
            else
                vec=obj.expr_saved;
            end
        end
        function expr=getType(obj)
            expr=obj.type;
        end
        function obj3=plus(obj1,obj2)
            assert((isa(obj1,'Evaluable') | isa(obj1,'double')) & (isa(obj2,'Evaluable') | isa(obj2,'double')),'Invalid use of PLUS - elements are not compatible (PET class: PrExpression)');
            if isa(obj1,'Evaluable') && isa(obj2,'Evaluable')
                obj3=Expression(obj1,1,obj2,1);
            elseif isa(obj1,'double')
                obj3=PrExpression(obj2.expression1,obj2.expression2,obj2.coef,obj2.pcoef+obj1);
            else %if isa(obj1,'double'))
                obj3=PrExpression(obj1.expression1,obj1.expression2,obj1.coef,obj1.pcoef+obj2);
            end
        end
        function obj3=minus(obj1,obj2)
            obj3=obj1+(-obj2);
        end
        function obj3=uminus(obj1)
            obj3=PrExpression(obj1.expression1,obj1.expression2,-obj1.coef,-obj1.pcoef);
        end
        function obj3=mtimes(a,b)
            assert((isa(a,'double') && isa(b,'PrExpression')) || (isa(b,'double') && isa(a,'PrExpression')),'Invalid use of MTIMES - elements are not compatible (PET class: PrExpression)');
            if isa(a,'double') && isa(b,'PrExpression')
                obj3=PrExpression(b.expression1,b.expression2,a*b.coef,a*b.pcoef);
            else
                obj3=PrExpression(a.expression1,a.expression2,b*a.coef,b*a.pcoef);
            end
        end
        function obj3=mrdivide(a,b)
            assert(isa(b,'double') && isa(a,'PrExpression'),'Invalid use of MRDIVIDE - elements are not compatible (PET class: PrExpression)');
            assert(b~=0,'Warning: division by 0');
            obj3=PrExpression(a.expression1,a.expression2,a.coef/b,a.pcoef/b);
        end
    end
end