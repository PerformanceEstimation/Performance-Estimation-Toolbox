function demo1
clc;
msg_length=90;
maketitle('Performance Estimation Toolbox (PESTO) -- EXAMPLE 1: a simple gradient method',msg_length,2)
fprintf('In this demo, we illustrate the use of PESTO for performing a worst-case \nanalysis of a gradient method (Examples/GradientMethod.m)\n')
waitfor;
fprintf('Consider the problem of minimizing a function F:\n\n min F(x),\n\nwhich is L-smooth and mu-strongly convex (say: L=1, mu=0.1).\n');
waitfor;
fprintf('We use the standard fixed-step gradient method\n\n x_{k+1}=x_k-h*g_k (g_k is the gradient of F at x_k),\n\nwith the fixed step size h=1/L.\n');
waitfor;
fprintf('The next lines show how to compute the worst-case objective function accuracy F(x_N)-F(xs) \n(xs is the optimum) assuming we start at some x_0, satisfying the condition\n\n||x_0-xs||^2<=1.\n');
waitfor;
maketitle('PESTO''s philosophy',msg_length,1)
waitfor;
fprintf('The goal is to provide the researchers an easy access to the performance estimation \nmethodology. The toolbox allows computing worst-case performances of simple first-order \nmethods by writing the algorithms (nearly) as if we were implementing them.\n');
waitfor;
maketitle('Initializing an empty performance estimation problem',msg_length,1)
waitfor;
disp('% (0) Initialize an empty PEP: P');
disp('>> P=pep();');
P=pep();
waitfor;
maketitle('Setting up the objective function F(x)',msg_length,1)
waitfor;
fprintf('%% (1) Setting up the objective function \n%%     We set the values of mu and L within a structure ''param'' :\n');
disp('>> param.mu=.1;   % Strong convexity parameter');
disp('>> param.L=1;     % Smoothness parameter');
param.mu=.1;
param.L=1;
waitfor;
fprintf('%% F is the L-smooth mu-strongly convex objective function:\n');
disp('>> F=P.DeclareFunction(''SmoothStronglyConvex'',param);');
waitfor;
fprintf('For more details about smooth strongly convex function, ''help SmoothStronglyConvex''.\nFor other functional classes, see User''s Manual, ''help pep'' or ''help pesto''.\n');
F=P.DeclareFunction('SmoothStronglyConvex',param);
waitfor;
maketitle('Setting up the starting point and initial condition',msg_length,1)
waitfor;
disp('% (2) Set up the starting point and initial condition');
disp('>> x0=P.StartingPoint(); % starting point x0');
waitfor;
disp('>> [xs,fs]=F.OptimalPoint(); % xs optimal and fs=F(xs)');
waitfor;
disp('>> P.InitialCondition((x0-xs)^2<=1); % ||x0-xs||^2<= 1');
% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1
waitfor;
maketitle('Setting up the algorithm',msg_length,1);
waitfor;
fprintf('As previously underlined, the idea is to write the algorithm is the same way as we would\nimplement it.\n\n');
waitfor;
waitfor
disp('% (3) Fixed-step gradient scheme; set step size and number of iterations')
disp('');
disp('>> h=1/param.L; % step size');
disp('>> N=10; % number of iterations');
h=1/param.L;
N=10;
waitfor;
disp('% Algorithm');
disp('% We use the routine ''oracle'' associated to F to evaluate function values and gradients');
disp('>> x=x0;');
disp('>> for i=1:N');
disp('>>     [g,f]=F.oracle(x); % g=grad F(x), f=F(x)');
disp('>>     x=x-h*g;');
disp('>> end');
disp('>> xN=x; %xN is the last iterate');
x=x0;
for i=1:N
    [g,f]=F.oracle(x);		% g=grad F(x), f=F(x)
    x=x-h*g;
end
xN=x;
waitfor;
fprintf('Note that only simple arithmetic operations are used for generating the algorithm:\n')
fprintf('(a) sums and differences of elements of the same natures:\n   - elements of the decision space on the one hand: x''s and g''s, and\n   - scalar values on the other hand (e.g.,f''s);\nand (b) scalar products between elements of the decision space: e.g., x*x, x*g, g*g.\n(we focus on the so called ''linearly Gram-representable'' cases)');
waitfor;
maketitle('Setting up the performance measure',msg_length,1);
waitfor
disp('% (4) Set up the performance measure');
disp('>> [g,f]=F.oracle(xN); % g=grad F(xN), f=F(xN)');
disp('>> P.PerformanceMetric(f-fs); % Worst-case evaluated as F(xN)-F(xs)');

% (4) Set up the performance measure
[g,f]=F.oracle(xN);                % g=grad F(xN), f=F(xN)
P.PerformanceMetric(f-fs); % Worst-case evaluated as F(xN)-F(xs)

waitfor;
maketitle('Solve the PEP',msg_length,1);
waitfor;
disp('% (5) Solve the PEP');
disp('>> P.solve()');
fprintf('\n');
maketitle('Solving the PEP: computations',msg_length,1);
waitfor;

% (5) Solve the PEP
P.solve()
waitfor;
maketitle('Exploiting the outputs',msg_length,1);
waitfor;
fprintf('In order to use the outputed values, the main idea is to use the command ''double()''.\nThis command allows to evaluate the values of the different elements of the worst-case\nidentified by the solver (e.g., function values, gradients, coordinates,...).\n\n');
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
fprintf('Note that the performance estimation problem is solved in terms of a Gram matrix\ncontaining the scalar products between x''s and g''s. For obtaining the x''s and g''s, \nwe use an approximate Cholesky decomposition of the matrix (by discarding \nthe eigenvalues smaller than 1e-9).\n')
waitfor;
fprintf('End of DEMO 1: for further details, we refer to demo2 and to the User''s Manual.\n');
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
