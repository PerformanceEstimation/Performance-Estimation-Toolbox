function demo2
clear all; clc;
msg_length=90; Nmax=10;
maketitle('Performance Estimation Toolbox (PESTO) -- EXAMPLE 2: comparing subgradient methods',msg_length,2)
fprintf('In this demo, we illustrate the use of PESTO for performing a worst-case\nanalysis of a subgradient method (Examples/SubgradientMethod.m).\n')
waitfor;
fprintf('After that, we show how to use this worst-case computation problem to compare different\nsimple standard step size policies (for N=1,...,%d).\n',Nmax);
waitfor;
fprintf('\nConsider the problem of minimizing a function F:\n\n min F(x),\n\nwhich is convex with bounded subgradients (Lipschitz): ||g||<=L (say, with L=1).\n');
waitfor;
fprintf('We use a subgradient method with step size policy h \n\n x_{k+1}=x_k-h_k g_k (g_k is a subgradient of F at x_k).\n\n');
waitfor;
fprintf('The next lines show how to compute the worst-case objective function accuracy of the\nbest iterate min_{0<=i<=N} (F(x_i)-F(xs)) (xs is the optimum) assuming we start at some\nx_0, satisfying the condition\n\n ||x_0-xs||^2<=1.\n');
waitfor;
maketitle('Embedding the PEP within a function',msg_length,1)
waitfor;
fprintf('For comparing the different policies, we embed the PEP within the following function.\n\n');
disp('>> function WCPerformance=PEPforSubgradient(h)');
disp('>>  ... % PEP code here (described hereafter)');
disp('>> end');
fprintf('\nwhere h is a vector of size N (number of iterations) containing all the step sizes.\n');
waitfor;
maketitle('Initializing an empty performance estimation problem',msg_length,1)
waitfor;
disp('% (0) Initialize an empty PEP');
disp('>> P=pep();');
waitfor;
maketitle('Setting up the objective function F(x)',msg_length,1)
waitfor;
fprintf('%% (1) Set up the objective function \n%%     We set the Lipschitz constant within a structure ''param'' :\n');
disp('>> param.R=1;	% ''radius''-type constraint on the subgradient norms: ||g||<=1');
waitfor;
fprintf('%% F is the objective function:\n');
disp('>> F=P.DeclareFunction(''SmoothConvexBoundedGradient'',param);');
waitfor;
fprintf('Details on convex functions with bounded subgradients, ''help SmoothConvexBoundedGradient''.\nFor other functional classes, see User''s Manual, ''help pep'' or ''help pesto''.\n');
waitfor;
maketitle('Setting up the starting point and initial condition',msg_length,1)
waitfor;
disp('% (2) Set up the starting point and initial condition');
disp('>> x0=P.StartingPoint();	% starting point x0');
waitfor;
disp('>> [xs,fs]=F.OptimalPoint();	% xs optimal and fs=F(xs)');
waitfor;
disp('>> P.InitialCondition((x0-xs)^2<=1); %  ||x0-xs||^2<= 1');

waitfor;
maketitle('Setting up the algorithm and performance measure',msg_length,1);
waitfor;
fprintf('As in the previous example, the idea is to write the algorithm is the same way as we would\nimplement it.\n\n');
waitfor;
disp('% (3) Algorithm and (4) performance measure');
disp('>> N=length(h); % number of iterations (remember: h is the input containing all step sizes)');
fprintf('\n');
disp('% Note: the worst-case performance measure used in the PEP is the');
disp('%       min_i (PerformanceMetric_i) (i.e., the best value among all');
disp('%       performance metrics added into the problem). Here, we use it');
disp('%       in order to find the worst-case value for min_i [F(x_i)-F(xs)].');
waitfor;
fprintf('\n');
disp('>> x=x0;');
disp('>> for i=1:N');
disp('>>     [g,f]=F.oracle(x);		% g=grad F(x), f=F(x)');
disp('>>     P.PerformanceMetric(f-fs);');
disp('>>     x=x-h(i)*g;');
disp('>> end');
disp('>> xN=x;');
disp('>> ');
disp('>> [~,f]=F.oracle(xN);');
disp('>> P.PerformanceMetric(f-fs);');
waitfor;

maketitle('Solve the PEP',msg_length,1);
waitfor;
disp('% (5) Solve the PEP');
disp('>> verbose=0;');
disp('>> sol=P.solve(verbose); % solves the PEP and disables the verbose mode.');
disp('>> WCPerformance=sol.WCperformance;');
disp('');

waitfor;

maketitle('Comparison of the step size policies',msg_length,1);
waitfor;
fprintf('Evaluating the worst-case of three policies for N=1,...,%d; this may take a few minutes.\n',Nmax);
fprintf('\nPolicy 1: h_i=1/sqrt(i+1) (green). ')
perf_policy1=zeros(Nmax,1);
msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.\n', 0,Nmax);
fprintf(msg);
prevlength = numel(msg);
for i=1:Nmax
    msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.\n', i,Nmax);
    fprintf(repmat('\b', 1, prevlength));
    fprintf(msg);
    h=1./sqrt((1:i)+1);
    prevlength = numel(msg);
    perf_policy1(i)=PEPforSubgradient(h);
end
scrsz = get(0,'ScreenSize');
 figure('Position',[1+scrsz(3)/2 1+scrsz(3)/2 scrsz(3)/2 scrsz(4)/2])
hold on; xlabel('Iteration count'); ylabel('min_i (F(x_i)-F(x_*))');
plot(1:Nmax,perf_policy1,'-g','linewidth',2);
fprintf('\nPolicy 2: h_i=1/(i+1) (blue). ')
perf_policy2=zeros(Nmax,1);
msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.', 0,Nmax);
fprintf(msg);
prevlength = numel(msg);
for i=1:Nmax
    msg = sprintf('Worst-case computation for N=%d out of Nmax=%d done.\n', i,Nmax);
    fprintf(repmat('\b', 1, prevlength));
    fprintf(msg);
    h=1./((1:i)+1);
    prevlength = numel(msg);
    perf_policy2(i)=PEPforSubgradient(h);
end
plot(1:Nmax,perf_policy2,'-b','linewidth',2);
fprintf('\nPolicy 3: h_i=1/sqrt(N+1) (constant, optimal step sizes) (red). ')
perf_policy3=zeros(Nmax,1);
msg = sprintf('Worst-case computation \nfor N=%d out of Nmax=%d done.\n', 0,Nmax);
fprintf(msg);
prevlength = numel(msg);
for i=1:Nmax
    msg = sprintf('Worst-case computation \nfor N=%d out of Nmax=%d done.\n', i,Nmax);
    fprintf(repmat('\b', 1, prevlength));
    fprintf(msg);
    h=ones(i,1)*1/sqrt(i+1);
    prevlength = numel(msg);
    perf_policy3(i)=PEPforSubgradient(h);
end
plot(1:Nmax,perf_policy3,'--r','linewidth',2);

waitfor;
fprintf('End of DEMO 2: for further details, we refer to demo3 and to the User''s Manual.\n');
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
function WCPerformance=PEPforSubgradient(h)
% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.R=1;     % Lipschitz constant (radius-style)

F=P.DeclareFunction('SmoothConvexBoundedGradient',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1);

% (3) Algorithm and (4) performance measure
N=length(h); % number of iterations

% Note: the worst-case performance measure used in the PEP is the 
%       min_i (PerformanceMetric_i) (i.e., the best value among all
%       performance metrics added into the problem). Here, we use it
%       in order to find the worst-case value for min_i [F(x_i)-F(xs)]
x=x0;
for i=1:N
    [g,f]=F.oracle(x);
    P.PerformanceMetric(f-fs);
    x=x-h(i)*g;
end
xN=x;

[~,f]=F.oracle(xN);
P.PerformanceMetric(f-fs);

% (5) Solve the PEP
sol=P.solve(0);
WCPerformance=sol.WCperformance;
end