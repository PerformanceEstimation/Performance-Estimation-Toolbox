function demo2
clear all; clc;
msg_length=90;
maketitle('Performance Estimation Toolbox (PET) -- DEMO 2 (Subgradient Method)',msg_length,2)
fprintf('In this demo, we illustrate the use of the PET toolbox for performing\na worst-case analysis of a subgradient method (Examples/SubgradientMethod.m)')
waitfor;
fprintf('Consider the problem of minimizing a function F:\n\n min F(x),\n\nwhich is convex with bounded subgradients (Lipschitz): ||g||<=L (with L=1).');
waitfor;
fprintf('We use a subgradient method with step size policy gamma_k \n\n x_{k+1}=x_k-gamma_k g_k (g_k is a subgradient of F at x_k),\n\n');
waitfor;
fprintf('The next lines show how to compute the worst-case objective function accuracy of the\nbest iterate min_{0<=i<=N} (F(x_i)-F(xs)) (xs is the optimum) assuming we start at some\nx_0, satisfying the condition\n\n ||x_0-xs||^2<=1.\n');
waitfor;

Nmax=10;
fprintf('After that, we show how to use this worst-case computation problem to compare different\nsimple standard step size policies (for N=1,...,%d).\n',Nmax)
waitfor;

maketitle('Embedding the PEP within a function',msg_length,1)
waitfor;
fprintf('For comparing the different policies, we embed the PEP within the following function.\n\n');
disp('>> function WCPerformance=PEPforSubgradient(gamma)');
disp('>>  ... % PEP code here (described hereafter)');
disp('>> end');
fprintf('\nwhere gamma is a vector containing all the step sizes\n');
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
disp('>> out=my_pep.solve();');
disp('>> WCPerformance=sol.WCperformance;');
disp('');

waitfor;

maketitle('Comparison of the step size policies',msg_length,1);
waitfor;
fprintf('Evaluating the worst-case of three policies for N=1,...,%d; this may take a few minutes.\n',Nmax);
fprintf('\nPolicy 1: gamma_i=1/sqrt(i+1) (green). ')
perf_policy1=zeros(Nmax,1);
msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.\n', 0,Nmax);
fprintf(msg);
prevlength = numel(msg);
for i=1:Nmax
    msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.\n', i,Nmax);
    fprintf(repmat('\b', 1, prevlength));
    fprintf(msg);
    gamma=1./sqrt((1:i)+1);
    prevlength = numel(msg);
    perf_policy1(i)=PEPforSubgradient(gamma);
end
scrsz = get(0,'ScreenSize');
 figure('Position',[1+scrsz(3)/2 1+scrsz(3)/2 scrsz(3)/2 scrsz(4)/2])
hold on; xlabel('Iteration count'); ylabel('min_i (F(x_i)-F(x_*))');
plot(1:Nmax,perf_policy1,'-g','linewidth',2);
fprintf('\nPolicy 2: gamma_i=1/(i+1) (blue). ')
perf_policy2=zeros(Nmax,1);
msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.', 0,Nmax);
fprintf(msg);
prevlength = numel(msg);
for i=1:Nmax
    msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.\n', i,Nmax);
    fprintf(repmat('\b', 1, prevlength));
    fprintf(msg);
    gamma=1./((1:i)+1);
    prevlength = numel(msg);
    perf_policy2(i)=PEPforSubgradient(gamma);
end
plot(1:Nmax,perf_policy2,'-b','linewidth',2);
fprintf('\nPolicy 3: gamma_i=1/sqrt(N+1) (constant, optimal step sizes) (red). ')
perf_policy3=zeros(Nmax,1);
msg = sprintf('Worst-case computation \nfor N=%d out of Nmax=%d done.\n', 0,Nmax);
fprintf(msg);
prevlength = numel(msg);
for i=1:Nmax
    msg = sprintf('Worst-case computation \nfor N=%d out of Nmax=%d done.\n', i,Nmax);
    fprintf(repmat('\b', 1, prevlength));
    fprintf(msg);
    gamma=ones(i,1)*1/sqrt(i+1);
    prevlength = numel(msg);
    perf_policy3(i)=PEPforSubgradient(gamma);
end
plot(1:Nmax,perf_policy3,'--r','linewidth',2);

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
function WCPerformance=PEPforSubgradient(gamma)
% (0) Initialize an empty PEP
my_pep=pet();

% (1) Set up the objective function
param.M1=Inf;	%
param.M2=1;     %
param.L=Inf;    % Smoothness parameter

F=my_pep.AddComponentObjective('SmoothConvexBoundedGradient',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=my_pep.GenStartingPoint();		 % x0 is some starting point
[xs,fs]=F.GetOptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)

my_pep.AddConstraint((x0-xs)*(x0-xs)<=1);

% (3) Algorithm
N=length(gamma);		% number of iterations
% gam=1/sqrt(N+1);		% step size

x=cell(N+1,1);
x{1}=x0;
for i=1:N
    [g,f]=F.oracle(x{i});
    my_pep.AddPerformanceConstraint(f-fs);
    x{i+1}=x{i}-gamma(i)*g;
end

[~,f]=F.oracle(x{N+1});
my_pep.AddPerformanceConstraint(f-fs);

% (5) Solve the PEP
sol=my_pep.solve();
WCPerformance=sol.WCperformance;
end