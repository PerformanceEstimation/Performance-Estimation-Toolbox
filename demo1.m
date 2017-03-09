function demo1
clear all; clc;
msg_length=90;
maketitle('Performance Estimation Toolbox (PET) -- DEMO 1 (Gradient Method)',msg_length,2)
fprintf('In this demo, we illustrate the use of the PET toolbox for performing\na worst-case analysis of a gradient method (Examples/GradientMethod.m)')
waitfor;
fprintf('Consider the problem of minimizing a function F:\n\n min F(x),\n\nwhich is L-smooth and mu-strongly convex (L=1, mu=0.1).');
waitfor;
fprintf('We use the standard fixed-step gradient method\n\n x_{k+1}=x_k-h*g_{k} (g_k is the gradient of F at x_k),\n\nwith the fixed step size gamma=1/L.\n');
waitfor;
fprintf('The next lines show how to compute the worst-case objective function accuracy F(x_N)-F(xs) \n(xs is the optimum) assuming we start at some x_0, satisfying the condition\n\n||x_0-xs||^2<=1.\n');
waitfor;
maketitle('Initializing an empty performance estimation problem',msg_length,1)
waitfor;
disp('% (0) Initialize an empty PEP');
disp('>> my_pep=pet();');
my_pep=pet();
waitfor;
maketitle('Setting up the objective function',msg_length,1)
waitfor;
fprintf('%% (1) Set up the objective function \n%%     We set the values of mu and L within a structure ''param'' :\n');
disp('>> param.mu=.1;   % Strong convexity parameter');
disp('>> param.L=1;     % Smoothness parameter');
param.mu=0;
param.L=1;
waitfor;
fprintf('%% F is the L-smooth mu-strongly convex objective function:\n');
disp('>> F=my_pep.AddComponentObjective(''SmoothStronglyConvex'',param);');
F=my_pep.AddComponentObjective('SmoothStronglyConvex',param);
waitfor;
maketitle('Setting up the starting point and initial condition',msg_length,1)
waitfor;
disp('% (2) Set up the starting point and initial condition');
disp('>> x0=P.GenStartingPoint();	% starting point x0');
waitfor;
disp('>> [xs,fs]=F.GetOptimalPoint();	% xs optimal and fs=F(xs)');
waitfor;
disp('>> my_pep.AddInitialCondition((x0-xs)^2<=1); %  ||x0-xs||^2<= 1');
% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
my_pep.AddInitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1
waitfor;
maketitle('Setting up the algorithm',msg_length,1);
waitfor
disp('% (3) Fixed-step gradient scheme; set step size and number of iterations')
disp('');
disp('>> gamma=1/param.L; % step size');
disp('>> N=10;            % number of iterations');
gamma=1/param.L;
N=10;
waitfor;
disp('% Algorithm');
disp('>> x=x0;');
disp('>> for i=1:N');
disp('>>     [g,f]=F.oracle(x);		% g=grad F(x), f=F(x)');
disp('>>     x=x-gam/param.L*g;');
disp('>> end');

x=x0;
for i=1:N
    [g,f]=F.oracle(x);		% g=grad F(x), f=F(x)
    x=x-gamma*g;
end
waitfor;
fprintf('Note that only simple arithmetic operations are for generating the algorithm:\n')
fprintf('(a) sums and difference of elements of the same natures:\n   - elements of the decision space on the one hand: x''s and g''s, and\n   - scalar values on the other hand (e.g.,f''s);\nand (b) scalar products between elements of the decision space: e.g., x*x, x*g, g*g.\n');
waitfor;
maketitle('Setting up the performance measure',msg_length,1);
waitfor
disp('% (4) Set up the performance measure');
disp('>> [g,f]=F.oracle(x);                % g=grad F(x), f=F(x)');
disp('>> my_pep.AddPerformanceConstraint(f-fs); % Worst-case evaluated as F(x)-F(xs)');

% (4) Set up the performance measure
[g,f]=F.oracle(x);                % g=grad F(x), f=F(x)
my_pep.AddPerformanceConstraint(f-fs); % Worst-case evaluated as F(x)-F(xs)

waitfor;
maketitle('Solve the PEP',msg_length,1);
waitfor;
disp('% (5) Solve the PEP');
disp('>> my_pep.solve()');
disp('');
maketitle('Solving the PEP: computations',msg_length,1);
waitfor;

% (5) Solve the PEP
my_pep.solve(1)
waitfor;
maketitle('Exploiting the outputs',msg_length,1);
waitfor;
fprintf('In order to use the outputed value, the main idea is to use the command ''double()'',\n\n');
fprintf('which evaluates the value of the variable from the solution of the PEP.\n')
fprintf('Example 1: evaluate the value of the worst-case measure:\n');
disp('>> double(f-fs)');
double(f-fs)
waitfor;
fprintf('Example 2: evaluate the vector x-xs:\n');
disp('>> double(x-xs)');
double(x-xs)
waitfor;
fprintf('Example 3: evaluate the distance ||x-xs||^2:\n');
disp('>> double((x-xs)^2)');
double((x-xs)^2)
waitfor;
fprintf('Example 4: evaluate the gradient norm ||g||^2:\n');
disp('>> double(g^2)');
double(g^2)
waitfor;
fprintf('End of DEMO 1: for further details, we refer to demo2 and to the User Manual.\n');
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
