classdef (Abstract) functionHandler < handle
    properties (Access=private)
        
        expr_list_others;
        list_size_others;
        names_list_others;
    end
    methods
        function obj=functionHandler()
            obj.list_size_others=0;
            obj.expr_list_others=cell(0,1);
            obj.names_list_others=cell(0,1);
            functionHandler.CountActive(1);
        end
        function delete(obj)
            functionHandler.CountActive(-1);
        end
        function obj=AddConstraint(obj,expr,name)
            assert(isa(expr,'Constraint'),'Invalid constraint');
            if nargin < 3
                obj.list_size_others=obj.list_size_others+1;
                obj.expr_list_others{obj.list_size_others,1}=expr;
                obj.names_list_others{obj.list_size_others,1}='Unnamed constraint';
            elseif nargin == 3
                assert(ischar(name),'Invalid name in constraint (must be array of characters; use single quotes)');
                obj.list_size_others=obj.list_size_others+1;
                obj.expr_list_others{obj.list_size_others,1}=expr;
                obj.names_list_others{obj.list_size_others,1}=name;
            end
        end
        function g=evaluate(obj,expr,tag) % SAME AS GRADIENT
            if nargin>=3
                assert(isa(tag,'char'),'Value call: second argument must be a tag (string)');
                spec=tag;
            else
                spec='';
            end
            [g,~]=obj.oracle(expr,spec);
        end
        function g=gradient(obj,expr,tag)
            if nargin>=3
                assert(isa(tag,'char'),'Value call: second argument must be a tag (string)');
                spec=tag;
            else
                spec='';
            end
            [g,~]=obj.oracle(expr,spec);
        end
        function g=subgradient(obj,expr,tag)
            if nargin>=3
                assert(isa(tag,'char'),'Value call: second argument must be a tag (string)');
                spec=tag;
            else
                spec='';
            end
            [g,~]=obj.oracle(expr,spec);
        end
        function f=value(obj,expr,tag)
            if nargin>=3
                assert(isa(tag,'char'),'Value call: second argument must be a tag (string)');
                spec=tag;
            else
                spec='';
            end
            [~,f]=obj.oracle(expr,spec);
        end
        function [cons,names]=collect(obj)
            cons=[];
            for i=1:obj.list_size_others
                lexpr=obj.expr_list_others{i,1}.Eval();
                cons=cons+lexpr;
            end
            names=obj.names_list_others;
        end
        function obj3=plus(obj1,obj2)
            assert(isa(obj1,'functionHandler') && isa(obj2,'functionHandler'));
            %             obj3=CompositeFunction(obj1,obj2);
            [table1]    = DecomposeCompositeFunction(obj1, 1);
            [table2]    = DecomposeCompositeFunction(obj2, 1);
            [tableTotal] = SummarizedCompositeFunction([table1; table2]);
            obj3 = FormulateCompositeFunctions(tableTotal);
        end
        function obj3=mtimes(obj1,obj2)
            assert((isa(obj1,'functionHandler') && isa(obj2,'double')) || (isa(obj1,'double') && isa(obj2,'functionHandler')),'Invalid multiplication');
            if isa(obj2,'double')
                obj3=Scaled_Function(obj1,obj2);
            else
                obj3=Scaled_Function(obj2,obj1);
            end
        end
        function obj3=mrdivide(obj1,obj2)
            assert(isa(obj1,'functionHandler') && isa(obj2,'double'),'Invalid division');
            obj3=Scaled_Function(obj1,1/obj2);
        end
        function obj3=minus(obj1,obj2)
            assert(isa(obj1,'functionHandler') && isa(obj2,'functionHandler'));
            obj3=obj1 + ((-1)*obj2);
        end
        function obj3=uminus(obj1)
            assert(isa(obj1,'functionHandler'));
            obj3= (-1)*obj1;
        end
    end
    methods (Abstract)
        OptimalPoint(obj,tag)
        oracle(obj,x,tag)
    end
    methods (Static)
        function out=CountActive(add)
            persistent nbEval;
            if isempty(nbEval)
                nbEval=0;
            end
            if nbEval+add>=0
                nbEval=nbEval+add;
            end
            out=nbEval;
        end
    end
end