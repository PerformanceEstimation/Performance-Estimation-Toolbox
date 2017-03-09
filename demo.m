function demo
clear all; clc;
msg_length=90;
maketitle('Performance Estimation Toolbox (PET) -- DEMONSTRATION MODULES',msg_length,2)
fprintf('Three demonstration modules are provided within the toolbox:\n\n')
fprintf('(1) demo1 provides the most basic example on how to use the toolbox for evaluating the\nworst-case performance of a gradient method for smooth strongly convex minimization.\n\n');
fprintf('(2) demo2 provides an example on how to use the toolbox for comparing the worst-case\nperformances of different step size policies for subgradient methods in non-smooth \nconvex minimization.\n\n');
fprintf('(3) demo3 provides an example on how to use the toolbox for comparing the worst-case\nperformances of FISTA with the best corresponding theoretical guarantee.\n\n');
end
function maketitle(msg,width,nblines)
if nargin<3
    nblines=1;
end
fprintf('\n');
for i=1:nblines
    fprintf([repmat('*',1,width) '\n']);
end
msg_l=numel(msg);
left_margin=floor((width-msg_l)/2);
right_margin=width-msg_l-left_margin;
for i=1:nblines-1
    fprintf(['*' repmat(' ',1,width-2) '*\n']);
end
msg_mod=['*' repmat(' ',1,left_margin-1) msg repmat(' ',1,right_margin-1) '*\n'];
fprintf(msg_mod)
for i=1:nblines-1
    fprintf(['*' repmat(' ',1,width-2) '*\n']);
end
for i=1:nblines
    fprintf([repmat('*',1,width) '\n']);
end
fprintf('\n');
end
