classdef Point < Evaluable
    properties (GetAccess=private)
        number;
        expression;
        expr_saved;
        when_saved;
    end
    methods
        function obj=Point(type,expr)
            assert((strcmp(type,'Point') | strcmp(type,'Function value')),'Type must be Point or Function value (PEsTo class: Point)');
            obj.type=type;
            obj.expr_saved=[];
            obj.when_saved=0;
            if nargin >= 2
                assert((isa(expr,'Evaluable') | isa(expr,'double')) & (strcmp(type,'Point') | strcmp(type,'Function value')),'Assignment is not valid (PEsTo class: Point)');
                if isa(expr,'Evaluable')
                    obj.expression=expr;
                else
                    assert(expr==0,'Assignment is not valid (PEsTo class: Point)');
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
                        F=Evaluable.SetGetFunc();
                        vec=vec.'*F;
                    end
                elseif obj.expression==0
                    dim=1;
                    if strcmp(obj.type,'Point')
                        dim=Point.SetGetPts(obj.type,0);
                    end
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
        function value=double(obj)
            sol=Evaluable.Solved();
            assert(sol.status==1,'This PEP has not been solved yet');
            dim_when_solved=length(sol.P);
            assert(dim_when_solved==Point.GetSize('Point'),'This PEP was modified after being solved: solve the new PEP before evaluating its solution');
            dim_when_solved=length(sol.F);
            assert(dim_when_solved==Point.GetSize('Function value'),'This PEP was modified after being solved: solve the new PEP before evaluating its solution');
            if obj.when_saved==0
                obj.Eval();
            end
            if strcmp(obj.type,'Function value')
                value=double(obj.expr_saved);
            else
                value=sol.P*obj.expr_saved;
            end
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
            assert(strcmp(type,'Point') | strcmp(type,'Function value'),'Type must be Point or Function value (PEsTo class: Point)');
            out=Point.SetGetPts(type,0);
        end
        function out=Reset(time)
            persistent reset_time;
            if nargin==1
                type='Point';
                out=Point.SetGetPts(type,0);
                Point.SetGetPts(type,-out);
                type='Function value';
                out=Point.SetGetPts(type,0);
                Point.SetGetPts(type,-out);
                reset_time=time;
            end
            out=reset_time;
            
        end
    end
end
