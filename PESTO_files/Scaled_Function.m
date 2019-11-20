classdef Scaled_Function < functionHandler
    properties (GetAccess=private)
        func;
        scaling; % this function is equal to func * scaling
    end
    methods
        function obj=Scaled_Function(f1,scaling)
            assert(isa(f1,'functionHandler') && isa(scaling,'double'),'Invalid statement');
            obj.func=f1;
            obj.scaling=scaling;
        end
        function obj=AddComponent(obj,x,g,f,spec)
            assert(strcmp(x.getType(),'Point') & strcmp(g.getType(),'Point') & strcmp(f.getType(),'Function value'),'Wrong type representation');
            g1=g/obj.scaling;
            f1=f/obj.scaling;
            if nargin < 5
                spec = '';
            end
            obj.func.AddComponent(x,g1,f1,spec);
        end
        function obj=AddConstraint(obj,expr)
            assert(isa(expr,'Constraint'),'Invalid initial condition');
            obj.func.AddConstraint(expr);
        end
        function cons=GetInterp(obj)
            cons=[];
            cons=cons+obj.funcGetInterp();
        end
        function [x,f]=OptimalPoint(obj,tag)
            if nargin < 2
                tag='';
            end
            [x,f1]=obj.func.OptimalPoint(tag);
            f=f1*obj.scaling;
        end
        function [x,f]=GetOptimalPoint(obj,tag)
            fprintf('GetOptimalPoint is deprecated, consider using OptimalPoint instead\n');
            if nargin ==2
                [x,f]=obj.OptimalPoint(tag);
            else
                [x,f]=obj.OptimalPoint();
            end
        end
        function [g, f]=oracle(obj,x,tag)
            assert(isa(x,'char') | isa(x,'Evaluable'),'Oracle call: x must either be a tag or a point');
            if nargin>=3
                assert(isa(tag,'char'),'Oracle call: second argument must be a tag (string)');
                [g1, f1]=obj.func.oracle(x, tag);
            else
                [g1, f1]=obj.func.oracle(x);
            end
            g=g1*obj.scaling;
            f=f1*obj.scaling;
        end
        function [f, scaling] = getFunctions(obj)
            f = obj.func;
            scaling = obj.scaling; 
        end
        function disp(obj)
            fprintf('Scaled function\n');
        end
    end
end