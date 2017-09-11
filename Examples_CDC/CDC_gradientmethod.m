function wc=CDC_gradientmethod(N,gamma)

% Initialize an empty PEP:
P=pep();

param.mu=0.1;	% Strong convexity parameter
param.L=1;  % Smoothness parameter

% F is the objective function
F=P.DeclareFunction('SmoothStronglyConvex',param); 

% x0 is some starting point:
x0=P.StartingPoint();

% xs is an optimal point, and fs=F(xs):
[xs,fs]=F.OptimalPoint();

% Initial condition F(x0)-F(xs)<=1:
P.InitialCondition(F.value(x0)-fs<=1); 

x=x0;
for i=1:N
    % Overwrite the previous value of x:
    x=x-gamma(i)*F.gradient(x);
end
xN=x;

% Evaluate the gradient at the last iterate:
gN=F.gradient(xN);
% Worst-case evaluated as ||gN||^2:        
P.PerformanceMetric(gN^2); 

verbose=1;
out=P.solve(verbose);
wc=out.WCperformance;%double(gN^2)

end