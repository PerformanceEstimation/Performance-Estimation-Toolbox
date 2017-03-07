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
            assert(isa(expr,'Evaluable'),'Invalid initial condition');
            assert(strcmp(expr.getType(),'Point') || strcmp(expr.getType(),'Function value'),'Invalid initial condition');
            obj.list_size_others=obj.list_size_others+1;
            obj.expr_list_others{obj.list_size_others,1}=expr;
        end
        function cons=collect(obj)
            cons=[];
            for i=1:obj.list_size_others
                lexpr=ExpressionWrapper(obj.expr_list_others{i,1}.Eval(),obj.expr_list_others{i,1}.getType());
                cons=cons+(lexpr<=0);
            end
        end
    end
end