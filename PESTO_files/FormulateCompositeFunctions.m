function [F] = FormulateCompositeFunctions(table)


nb_elem_f   = size(table,1);
max_weight  = table{1,2};
max_ind     = 1;
for i = 2:nb_elem_f
    max_weight = max([max_weight abs(table{i,2})]);
    if max_weight == abs(table{i,2}), max_ind = i; end
end
F = table{max_ind,2} * table{max_ind,1};
for i = 1:nb_elem_f
    if i~=max_ind && table{i,2}~=0
        F = CompositeFunction(F,table{i,2} * table{i,1});
    end
end
end