function [simplified] = DecomposeCompositeFunction(compositeFunction, final)


% decompose the composite function into basic function, and rewrite
% the composite function so that every function in the sum appears only
% once.

current    = compositeFunction;
leaves     = cell(0,2);
if isa(current,'Scaled_Function')
    [f1, sc1] = getFunctions(current);
    [branch1] = DecomposeCompositeFunction(f1,0);
    size_branch = size(branch1, 1);
    for i = 1:size_branch
        leaves{i,1} = branch1{i,1};
        leaves{i,2} = branch1{i,2} * sc1;
    end
elseif isa(current,'functionClass')
    leaves{1,1} = current;
    leaves{1,2} = 1;
elseif isa(current,'CompositeFunction')
    [f1,f2]=getFunctions(current);
    [branch1] = DecomposeCompositeFunction(f1,0);
    [branch2] = DecomposeCompositeFunction(f2,0);
    leaves = [branch1; branch2];
else
    assert(false, 'Invalid call to DecomposeCompositeFunction');
end

if final
    [simplified] = SummarizedCompositeFunction(leaves);
else
    simplified = leaves;
end
end




