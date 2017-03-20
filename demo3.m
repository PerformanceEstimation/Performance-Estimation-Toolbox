function demo3
clear all; clc;
msg_length=90;
maketitle('Performance Estimation Toolbox (PEsTo) -- DEMO 3 (FISTA)',msg_length,2)
fprintf('In this demo, we illustrate the use of PEsTo for performing\na worst-case analysis of a subgradient method (Examples/SubgradientMethod.m)')
waitfor;
fprintf('Consider the problem of minimizing a function F:\n\n min F(x),\n\nwhich is convex with bounded subgradients (Lipschitz): ||g||<=L (with L=1).');
waitfor;
fprintf('We use a subgradient method with step size policy gamma_k \n\n x_{k+1}=x_k-gamma_k g_k (g_k is a subgradient of F at x_k),\n\n');
waitfor;
fprintf('The next lines show how to compute the worst-case objective function accuracy of the\nbest iterate min_{0<=i<=N} (F(x_i)-F(xs)) (xs is the optimum) assuming we start at some\nx_0, satisfying the condition\n\n ||x_0-xs||^2<=1.\n');
waitfor;
fprintf('After that, we show how to use this worst-case computation problem to compare different\nsimple standard step size policies.\n')
waitfor;
end
function waitfor
waitmsg=sprintf('\n\nPress any key to continue\n');
waitlength = numel(waitmsg);
disp(waitmsg)
pause
fprintf(repmat('\b', 1, waitlength));
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
