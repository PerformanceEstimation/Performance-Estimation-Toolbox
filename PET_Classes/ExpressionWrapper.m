classdef ExpressionWrapper < handle
    properties (GetAccess=private)
        vec;
        type; % 'Point', 'Function value', 'Scalar'
        dim;
    end
    methods
        function obj=ExpressionWrapper(vec,type)
            assert(size(vec,1)==size(vec,2) || size(vec,2)==1,'Invalid dimensions');
            assert(strcmp(type,'Point') | strcmp(type,'Function value') | strcmp(type,'Scalar'),'Wrong type representation');
            obj.vec=vec;
            obj.dim=size(vec);
            obj.type=type;
            if strcmp(type,'Function value')
                obj.dim=size(obj.vec);
                obj.type='Scalar';
            end
        end
        function expr=getType(obj)
            expr=obj.type;
        end
        function tf=isEqual(obj1,obj2)
            tf=min(obj1.Eval()-obj2.Eval()==0);
        end
        function expr=getDim(obj)
            expr=obj.dim;
        end
        function expr=Eval(obj)
            expr=obj.vec;
        end
        function obj3=plus(obj1,obj2)
            assert((isa(obj1,'ExpressionWrapper') && isa(obj2,'double')) || (isa(obj2,'ExpressionWrapper') && isa(obj1,'double')) || (isa(obj1,'ExpressionWrapper') && isa(obj2,'ExpressionWrapper')),'Invalid use of SUM - elements are not compatible (PET class: ExpressionWrapper)');
            if (isa(obj1,'ExpressionWrapper') && isa(obj2,'double'))
                assert(strcmp(obj1.getType(),'Scalar'),'Invalid use of SUM - elements are not of the same dimensions (PET class: ExpressionWrapper)');
                obj3=ExpressionWrapper(obj1.vec+obj2,obj1.type);
            elseif (isa(obj2,'ExpressionWrapper') && isa(obj1,'double'))
                assert(strcmp(obj2.getType(),'Scalar'),'Invalid use of SUM - elements are not of the same dimensions (PET class: ExpressionWrapper)');
                obj3=ExpressionWrapper(obj2.vec+obj1,obj2.type);
            else
                assert(obj1.dim(1)==obj2.dim(1) && obj1.dim(2)==obj2.dim(2) && strcmp(obj1.type, obj2.type),'Invalid use of SUM - elements are not of the same dimensions (PET class: ExpressionWrapper)')
                obj3=ExpressionWrapper(obj1.vec+obj2.vec,obj1.type);
            end
        end
        function obj3=minus(obj1,obj2)
            assert((isa(obj1,'ExpressionWrapper') && isa(obj2,'double')) || (isa(obj2,'ExpressionWrapper') && isa(obj1,'double')) || (isa(obj1,'ExpressionWrapper') && isa(obj2,'ExpressionWrapper')),'Invalid use of MINUS - elements are not compatible (PET class: ExpressionWrapper)');
            if (isa(obj1,'ExpressionWrapper') && isa(obj2,'double'))
                assert(strcmp(obj1.getType(),'Scalar'),'Invalid use of MINUS - elements are not of the same dimensions (PET class: ExpressionWrapper)');
                obj3=ExpressionWrapper(obj1.vec-obj2,obj1.type);
            elseif (isa(obj2,'ExpressionWrapper') && isa(obj1,'double'))
                assert(strcmp(obj2.getType(),'Scalar'),'Invalid use of MINUS - elements are not of the same dimensions (PET class: ExpressionWrapper)');
                obj3=ExpressionWrapper(obj2.vec-obj1,obj2.type);
            else
                assert(obj1.dim(1)==obj2.dim(1) && obj1.dim(2)==obj2.dim(2) && strcmp(obj1.type, obj2.type),'Invalid use of MINUS - elements are not of the same dimensions (PET class: ExpressionWrapper)')
                obj3=ExpressionWrapper(obj1.vec-obj2.vec,obj1.type);
            end
        end
        function obj3=uminus(obj1)
            obj3=ExpressionWrapper(-obj1.vec,obj1.type);
        end
        function obj3=mtimes(a,b)
            assert((isa(a,'ExpressionWrapper') && isa(b,'ExpressionWrapper')) || (isa(a,'double') && isa(b,'ExpressionWrapper')) || (isa(b,'double') && isa(a,'ExpressionWrapper')),'Invalid use of MTIMES - elements are not compatible (PET class: ExpressionWrapper)');
            if isa(a,'double') && isa(b,'ExpressionWrapper')
                obj3=ExpressionWrapper(a*b.vec, b.type);
            elseif isa(b,'double') && isa(a,'ExpressionWrapper')
                obj3=ExpressionWrapper(b*a.vec, a.type);
            else
                assert((isa(a,'ExpressionWrapper') && isa(b,'ExpressionWrapper')),'Invalid use of MTIMES - elements are not compatible (PET class: ExpressionWrapper)');
                assert(a.dim(1)==b.dim(1) && a.dim(2)==1 && b.dim(2)==1 && strcmp(a.type,'Point') && strcmp(b.type,'Point'),'Invalid use of MTIMES - elements do not have compatible dimensions (PET class: ExpressionWrapper)');
                G=ExpressionWrapper.SetGetGram();
                obj3=ExpressionWrapper(a.vec.'*G*b.vec,'Scalar');
            end
        end
        function obj3=mpower(a,b)
            assert(b==2,'Only squares are accepted in expressions (PET class: ExpressionWrapper)');
            obj3=mtimes(a,a);
        end
        function obj3=mrdivide(a,b)
            assert(isa(b,'double') && isa(a,'ExpressionWrapper'),'Invalid use of MRDIVIDE - elements are not compatible (PET class: ExpressionWrapper)');
            assert(b~=0,'Warning: division by 0');
            obj3=ExpressionWrapper(a.vec/b,a.type);
        end
        function expr=le(a,b)
            assert((isa(a,'double') && isa(b,'ExpressionWrapper')) || (isa(b,'double') && isa(a,'ExpressionWrapper')) || (isa(a,'ExpressionWrapper') && isa(b,'ExpressionWrapper')),'Invalid use of inequality constraint (PET class: ExpressionWrapper)')
            if (isa(a,'double') && isa(b,'ExpressionWrapper'))
                expr=(a-Eval(b)<=0);
            elseif (isa(b,'double') && isa(a,'ExpressionWrapper'))
                expr=(Eval(a)-b<=0);
            else
                expr=(Eval(a-b)<=0);
            end
        end
        function expr=ge(a,b)
            expr=le(b,a);
        end
        function expr=eq(a,b)
            assert((isa(a,'double') && isa(b,'ExpressionWrapper')) || (isa(b,'double') && isa(a,'ExpressionWrapper')) || (isa(a,'ExpressionWrapper') && isa(b,'ExpressionWrapper')),'Invalid use of equality constraint (PET class: ExpressionWrapper)')
            if (isa(a,'double') && isa(b,'ExpressionWrapper'))
                expr=(a-Eval(b)==0);
            elseif (isa(b,'double') && isa(a,'ExpressionWrapper'))
                expr=(Eval(a)-b==0);
            else
                expr=(Eval(a-b)==0);
            end
        end
    end
    methods (Static)

    end
end