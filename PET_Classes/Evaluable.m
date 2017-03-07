classdef (Abstract) Evaluable < handle
    methods
        function obj=Evaluable()
            Evaluable.update(1);
        end
        function disp(obj)
            if (strcmp(obj.getType(),'Point'))
                msg='PEP variable (type: decision variable or gradient)\n';
            else
                msg='PEP variable (type: function value)\n';
            end
            fprintf(msg);
        end
    end
    methods (Static)
        function out=update(upd)
            persistent last_update;
            if upd
                last_update=now;
            elseif isempty(last_update)
                last_update=0;
            end
            out=last_update;
        end
    end
end

