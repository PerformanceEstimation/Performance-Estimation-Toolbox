classdef Point < Evaluable
    properties (GetAccess=private)
        type;
        number;
        expression;
        expr_saved;
        when_saved;
    end
    methods
        function obj=Point(type,expr)
            assert((strcmp(type,'Point') | strcmp(type,'Function value')),'Type must be Point or Function value (PET class: Point)');
            obj.type=type;
            obj.expr_saved=[];
            obj.when_saved=0;
            if nargin >= 2
                assert((isa(expr,'Evaluable') | isa(expr,'double')) & (strcmp(type,'Point') | strcmp(type,'Function value')),'Assignment is not valid (PET class: Point)');
                if isa(expr,'Evaluable')
                    obj.expression=expr;
                else
                    assert(expr==0,'Assignment is not valid (PET class: Point)');
                    obj.expression=0;
                end
            else
                if (strcmp(type,'Point'))
                    obj.number=Point.SetGetPts(obj.type,1);
                elseif strcmp(type,'Function value')
                    obj.number=Point.SetGetPts(obj.type,1);
                end
                obj.expression=[];
            end
        end
        function vec=Eval(obj)
            if  isempty(obj.expr_saved) || (obj.when_saved<=Evaluable.update(0))
                if isempty(obj.expression)
                    dim=Point.SetGetPts(obj.type,0);
                    vec=zeros(dim,1);
                    vec(obj.number,1)=1;
                    if strcmp(obj.type,'Function value')
                        F=ExpressionWrapper.SetGetFunc();
                        vec=vec.'*F;
                    end
                elseif obj.expression==0
                    dim=Point.SetGetPts(obj.type,0);
                    vec=zeros(dim,1);
                else
                    vec=obj.expression.Eval();
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
            assert((isa(obj1,'Evaluable') | isa(obj1,'double')) & (isa(obj2,'Evaluable') | isa(obj2,'double')),'Invalid use of PLUS - elements are not compatible (PET class: Point)');
            if isa(obj2,'Evaluable') && isa(obj1,'Evaluable')
                obj1.getType()
                obj2.getType()
                assert(strcmp(obj1.getType(),obj2.getType()),'Invalid use of PLUS - elements are not compatible (PET class: Point)');
                obj3=Expression(obj1,1,obj2,1);
            elseif isa(obj1,'double')
                obj3=Expression(obj2,1,0,0,obj1);
            else %isa(obj1,'double')
                obj3=Expression(obj1,1,0,0,obj2);               
            end            
        end
        function obj3=minus(obj1,obj2)
            assert(strcmp(obj1.getType(),obj2.getType()),'Invalid use of MINUS - elements are not compatible (PET class: Point)');
            obj3=Expression(obj1,1,obj2,-1);
        end
        function obj3=uminus(obj1)
            obj3=Expression(obj1,-1,Point(obj1.getType(),0),0);
        end
        function obj3=mtimes(a,b)
            assert((isa(a,'double') || isa(a,'Point')) & (isa(b,'double') || isa(b,'Point')),'Invalid use of MTIMES - elements are not compatible (PET class: Point)');
            if isa(a,'double') && isa(b,'Point')
                obj3=Expression(b,a,Point(b.getType(),0),0);
            elseif isa(b,'double') && isa(a,'Point')
                obj3=Expression(a,b,Point(a.getType(),0),0);
            else
                obj3=PrExpression(a,b,1);
            end
        end
        function obj3=mpower(a,b)
            assert(b==2,'Only squares are accepted in expressions (PET class: Point)');
            obj3=mtimes(a,a);
        end
    end
    methods (Access=private, Static)
        function out=SetGetPts(type,add)
            assert(strcmp(type,'Point') | strcmp(type,'Function value'),'Wrong type representation');
            persistent nbPointsG;
            persistent nbPointsF;
            if isempty(nbPointsG)
                nbPointsG=0;
            end
            if isempty(nbPointsF)
                nbPointsF=0;
            end
            if (strcmp(type,'Point'))
                nbPointsG=nbPointsG+add;
                out=nbPointsG;
            else %(strcmp(type,'Function value'))
                nbPointsF=nbPointsF+add;
                out=nbPointsF;
            end
        end
    end
    methods (Static)
        function out=GetSize(type)
            assert(strcmp(type,'Point') | strcmp(type,'Function value'),'Type must be Point or Function value (PET class: Point)');
            out=Point.SetGetPts(type,0);
        end
    end
end
