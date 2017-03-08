classdef (Abstract) functionHandler < handle
    properties (Access=private)
        
        expr_list_others;
        list_size_others;
    end
    methods
        function obj=functionHandler()
            obj.list_size_others=0;
            obj.expr_list_others=cell(0,1);
        end
        function obj=AddConstraint(obj,expr)
            assert(isa(expr,'Constraint'),'Invalid initial condition');
            obj.list_size_others=obj.list_size_others+1;
            obj.expr_list_others{obj.list_size_others,1}=expr;
        end
        function cons=collect(obj)
            cons=[];
            for i=1:obj.list_size_others
                lexpr=obj.expr_list_others{i,1}.Eval();
                cons=cons+lexpr;
            end
        end
    end
end