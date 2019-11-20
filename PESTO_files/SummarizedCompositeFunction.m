function [simplified] = SummarizedCompositeFunction(leaves)

simplified = cell(0,2);
nb_elem    = size(leaves,1);
nb_elem_f  = 0;
for i = 1:nb_elem
    current_funct  = leaves{i,1};
    current_weight = leaves{i,2};
    if (leaves{i,2} ~= 0 || i == 1)
        nb_elem_f = nb_elem_f + 1;
        for j = i+1:nb_elem
            if leaves{j,1} == current_funct
                current_weight = current_weight + leaves{j,2};
                leaves{j,2}    = 0;
            end
        end
        simplified{nb_elem_f,1} = current_funct;
        simplified{nb_elem_f,2} = current_weight;
    end
end
end