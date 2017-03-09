classdef CompositeFunction < functionHandler
    properties (GetAccess=private)
        func1;
        func2;
    end
    methods
        function obj=CompositeFunction(f1,f2)
            assert(isa(f1,'functionHandler') && isa(f2,'functionHandler'),'Invalid statement');
            obj.func1=f1;
            obj.func2=f2;
        end
        function obj3=plus(obj1,obj2)
            assert(isa(obj1,'functionHandler') && isa(obj2,'functionHandler'));
            obj3=CompositeFunctions(obj1,obj2);
        end
        function obj=AddComponent(obj,x,g,f,spec)
            assert(strcmp(x.getType(),'Point') & strcmp(g.getType(),'Point') & strcmp(f.getType(),'Function value'),'Wrong type representation');
            g1=Point('Point');
            f1=Point('Function value');
            f2=f-f1;
            g2=g-g1;
            obj.func1.AddComponent(x,g1,f1,spec);
            obj.func2.AddComponent(x,g2,f2,spec);
        end
        function cons=GetInterp(obj)
            cons=[];
            cons=cons+obj.func1.GetInterp();
            cons=cons+obj.func2.GetInterp();
        end
        function [x,f]=GetOptimalPoint(obj,tag)
            x=Point('Point');
            g1=Point('Point');
            f1=Point('Function value');
            f2=Point('Function value');
            g2=-g1;
            if nargin < 2
                tag='';
            end
            obj.func1.AddComponent(x,g1,f1,tag);
            obj.func2.AddComponent(x,g2,f2,tag);
            f=f1+f2;
        end
        function [g, f]=oracle(obj,x,tag)
            assert(isa(x,'char') | isa(x,'Evaluable'),'Oracle call: x must either be a tag or a point');
            if nargin>=3
                assert(isa(tag,'char'),'Oracle call: second argument must be a tag (string)');
                spec=tag;
            else
                spec='';
            end
            if isa(x,'char')
                [g1, f1]=obj.func1.oracle(x);
                [g2, f2]=obj.func2.oracle(x);
            else
                g1=Point('Point');
                g2=Point('Point');
                f1=Point('Function value');
                f2=Point('Function value');
                obj.func1.AddComponent(x,g1,f1,spec);
                obj.func2.AddComponent(x,g2,f2,spec);
            end
            g=g1+g2;
            f=f1+f2;
        end
        function disp(obj)
            fprintf('Composite function\n');
        end
    end
end