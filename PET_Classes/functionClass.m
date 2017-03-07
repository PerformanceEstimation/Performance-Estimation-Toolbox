classdef functionClass < functionHandler
    properties (GetAccess=private)
        interp_list;
        list_size;
        interp_handler;
    end
    methods
        function obj=functionClass(InterpEval)
            obj.list_size=0;
            obj.interp_list=cell(0,1);
            assert(isa(InterpEval,'function_handle'),'Argument must be a handler');
            obj.interp_handler=InterpEval;
        end
        function obj=AddComponent(obj,x,g,f,spec)
            assert(strcmp(x.getType(),'Point') & strcmp(g.getType(),'Point') & strcmp(f.getType(),'Function value'),'Wrong type representation');
            obj.list_size=obj.list_size+1;
            if nargin > 4
                obj.interp_list{obj.list_size}.spec=spec;
            else
                obj.interp_list{obj.list_size}.spec='';
            end
            obj.interp_list{obj.list_size}.x=x;
            obj.interp_list{obj.list_size}.g=g;
            obj.interp_list{obj.list_size}.f=f;
        end
        function [x,f,g]=GetOptimalPoint(obj)
            x=Point('Point');
            g=Point('Point',0);
            f=Point('Function value');
            obj.AddComponent(x,g,f,'optimum');
        end
        function [g, f]=oracle(obj,x,reuse)
            if nargin > 2
                if reuse
                    for i=1:obj.list_size
                        if x==obj.interp_list{i}.x
                            g=obj.interp_list{i}.g;
                            f=obj.interp_list{i}.f;
                            found=1;
                            break;
                        end
                    end
                end
            else
                found=0;
            end
            if ~found
                g=Point('Point');
                f=Point('Function value');
                obj.AddComponent(x,g,f);
            end
            
        end
        function cons=GetInterp(obj)
            cons=[];
            for i=1:obj.list_size
                clear pt1;
                pt1.x=ExpressionWrapper(obj.interp_list{i}.x.Eval(),obj.interp_list{i}.x.getType());
                pt1.g=ExpressionWrapper(obj.interp_list{i}.g.Eval(),obj.interp_list{i}.g.getType());
                pt1.f=ExpressionWrapper(obj.interp_list{i}.f.Eval(),obj.interp_list{i}.f.getType());
                pt1.spec=obj.interp_list{i}.spec;
                for j=1:obj.list_size
                    clear pt2;
                    pt2.x=ExpressionWrapper(obj.interp_list{j}.x.Eval(),obj.interp_list{j}.x.getType());
                    pt2.g=ExpressionWrapper(obj.interp_list{j}.g.Eval(),obj.interp_list{j}.g.getType());
                    pt2.f=ExpressionWrapper(obj.interp_list{j}.f.Eval(),obj.interp_list{j}.f.getType());
                    pt2.spec=obj.interp_list{j}.spec;
                    new_cons=obj.interp_handler(pt1,pt2);
                    if ~isempty(new_cons)
                        cons=cons+new_cons;
                    end
                end
            end
        end
        function obj3=plus(obj1,obj2)
            assert(isa(obj1,'functionClass') && isa(obj2,'functionClass'));
            obj3=CompositeFunction(obj1,obj2);
        end
        function disp(obj)
            fprintf('Function, %d interpolation points\n',obj.list_size);
        end
    end
end