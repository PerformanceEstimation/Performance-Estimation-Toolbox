classdef Expression < Evaluable
    properties (GetAccess=private)
        expression1; % containing its dependence??
        coef1;
        expression2;
        coef2;
        pcoef;
        type;
        expr_saved;
        when_saved;
    end
    methods
        function obj=Expression(expr1,coef1,expr2,coef2,pcoef)
            assert(isa(expr1,'Evaluable') && isa(expr2,'Evaluable'),'Must be expression objects (PET class: Expression)');
            assert(strcmp(expr1.getType(),expr2.getType()),'Wrong type combination - incompatible for product expression (PET class: Expression)');
            
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
        function expr=getType(obj)
            expr=obj.type;
        end
        function obj3=plus(obj1,obj2)
            obj3=Expression(obj1,1,obj2,1);
        end
        function obj3=minus(obj1,obj2)
            obj3=Expression(obj1,1,obj2,-1);
        end
        function obj3=uminus(obj1)
            obj3=Expression(obj1.expression1,-obj1.coef1,obj1.expression2,-obj1.coef2);
        end
        function obj3=mtimes(a,b)
            assert((isa(a,'double') || isa(a,'Evaluable')) & (isa(b,'double') || isa(b,'Evaluable')),'Invalid use of MTIMES - elements are not compatible (PET class: Expression)');
            if isa(a,'double') && isa(b,'Expression')
                obj3=Expression(b.expression1,a*b.coef1,b.expression2,a*b.coef2,a*b.pcoef);
            elseif isa(b,'double') && isa(a,'Expression')
                obj3=Expression(a.expression1,b*a.coef1,a.expression2,b*a.coef2,b*a.pcoef);
            else
                obj3=PrExpression(a,b,1);
            end
        end
        function obj3=mpower(a,b)
            assert(b==2,'Only squares are accepted in expressions (PET class: Expression)');
            obj3=mtimes(a,a);
        end
        function obj3=mrdivide(a,b)
            assert(isa(b,'double') && isa(a,'Expression'),'Invalid use of MRDIVIDE - elements are not compatible (PET class: Expression)');
            assert(b~=0,'Warning: division by 0');
            obj3=Expression(a.expression1,a.coef1/b,a.expression2,a.coef2/b);
        end
    end
end