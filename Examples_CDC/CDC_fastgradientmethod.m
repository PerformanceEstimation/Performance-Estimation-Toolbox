function wc=CDC_fastgradientmethod(N,eps)
% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.L=1;      % Smoothness parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();             % x0 is some starting point
[xs,fs]=F.OptimalPoint();              % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
x=cell(N+1,1);% we store the iterates in a cell for convenience
x{1}=x0;
y=x0;
for i=1:N
    d  = inexactsubgradient(y,F,eps);
    x{i+1}=y-1/param.L*d;
    y=x{i+1}+(i-1)/(i+2)*(x{i+1}-x{i});
end

% (4) Set up the performance measure
[g,f]=F.oracle(x{N+1});                % g=grad F(x), f=F(x)
P.PerformanceMetric(f-fs); % Worst-case evaluated as F(x)-F(xs)

verbose=1;
out=P.solve(verbose);
wc=out.WCperformance;
end