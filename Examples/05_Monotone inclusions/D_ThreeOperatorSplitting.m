function D_ThreeOperatorSplitting
% In this example, we use a Three-Operator Splitting method for solving
% a monotone inclusion problem
%   find x st   0 \in Ax + Bx + Cx
% where A is maximally monotone, B is cocoercive and C is the gradient of 
% a smooth strongly convex function. We denote by JA and JB the respective
% resolvents of A and B.
%
% One iteration of the algorithm starting from point w is as follows:
%       x = JB( w )
%       y = JA( 2 * x - w - C x)
%       z = w - theta * ( x - y )
% and then we choose as the next iterate the value of z.
%
% Given two initial points w0 and w1, we show how to compute the worst-case
% contraction factor ||z0 - z1||/||w0 - w1|| obtained after doing one
% iteration of DRS from respectively w0 and w1. Note that we allow the 
% user to choose a stepsize alpha in the resolvent and for C.
%
% (0) Initialize an empty PEP
P=pep();

% (1) Set up the class of monotone inclusions
paramB.beta = 1;            % B is 1-cocoercive
paramC.L = 1; paramC.mu=.1; % C is the gradient of a 1-smooth .1-str convex

A = P.DeclareFunction('Monotone');
B = P.DeclareFunction('Cocoercive',paramB);
C = P.DeclareFunction('SmoothStronglyConvex',paramC);

% (2) Set up the starting points
w0=P.StartingPoint(); w1=P.StartingPoint();
P.InitialCondition((w0-w1)^2<=1);  % Normalize the initial distance ||w0-ws||^2 <= 1

% (3) Algorithm
alpha = .9;		% step size (in the resolvents)
theta = 1.3;    % overrelaxation

x0 = proximal_step(w0,B,alpha);
y0 = proximal_step(2*x0-w0-alpha*C.evaluate(x0),A,alpha);
z0 = w0-theta*(x0-y0);

x1 = proximal_step(w1,B,alpha);
y1 = proximal_step(2*x1-w1-alpha*C.evaluate(x1),A,alpha);
z1 = w1-theta*(x1-y1);

% (4) Set up the performance measure: ||z0-z1||^2
P.PerformanceMetric((z0-z1)^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double((z0-z1)^2)   % worst-case contraction factor
end