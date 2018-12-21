clear all; clc;
% In this example, we use a Douglas-Rachford splitting (DRS) 
% method for solving a monotone inclusion problem
%   find x st   0 \in Ax + Bx 
% where A is L-Lipschitz and monotone and B is (maximally) mu-strongly
% monotone. We denote by JA and JB the respective resolvents of A and B.
%
% One iteration of the algorithm starting from point w is as follows:
%       x = JB( w )
%       y = JA( 2 * x - w )
%       z = w - theta * ( x - y )
% and then we choose as the next iterate the value of z.
%
% Given two initial points w0 and w1, we show how to compute the worst-case
% contraction factor ||z0 - z1||/||w0 - w1|| obtained after doing one
% iteration of DRS from respectively w0 and w1.
% Note that we allow the user to choose a stepsize alpha in the resolvent.
%
% This setting is studied in
% (**) Ernest K. Ryu, Adrien B. Taylor, C. Bergeling, and P. Giselsson.
%      "Operator Splitting Performance Estimation: Tight contraction
%       factors and optimal parameter selection." arXiv:1812.00146, 2018.
%
% (0) Initialize an empty PEP
P=pep();

% (1) Set up the class of monotone inclusions
paramA.L  =  1; paramA.mu = 0; % A is 1-Lipschitz and 0-strongly monotone
paramB.mu = .1;                % B is .1-strongly monotone

A = P.DeclareFunction('LipschitzStronglyMonotone',paramA);
B = P.DeclareFunction('StronglyMonotone',paramB);

% (2) Set up the starting points
w0=P.StartingPoint();
w1=P.StartingPoint();
P.InitialCondition((w0-w1)^2<=1);  % Normalize the initial distance ||w0-ws||^2 <= 1

% (3) Algorithm
alpha = 1.3;		% step size (in the resolvents)
theta = .9;         % overrelaxation

x0 = proximal_step(w0,B,alpha);
y0 = proximal_step(2*x0-w0,A,alpha);
z0 = w0-theta*(x0-y0);

x1 = proximal_step(w1,B,alpha);
y1 = proximal_step(2*x1-w1,A,alpha);
z1 = w1-theta*(x1-y1);

% (4) Set up the performance measure: ||z0-z1||^2
P.PerformanceMetric((z0-z1)^2);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double((z0-z1)^2)   % worst-case contraction factor

% Results to be compared with WC below (from (**)):
L    = alpha*paramA.L; mu = alpha*paramB.mu;
C    = sqrt(((2*(theta-1)*mu+theta-2)^2+L^2*(theta-2*(mu+1))^2)/(L^2+1));
if theta*(theta+C)/(mu+1)^2/C * (C+mu*((2*(theta-1)*mu+theta-2)-L^2*...
        (theta-2*(mu+1)))/(L^2+1)) >=0
    WC   = ((theta+C)/2/(mu+1))^2;
elseif L<=1 && mu >= (L^2+1)/(L-1)^2 && theta<=-(2*(mu+1)*(L+1)*(mu+...
        (mu-1)*L^2-2*mu*L-1))/(mu+L*(L^2+L+1)+2*mu^2*(L-1)+mu*L*(1-(L-3)*L)+1)
    WC   = (1-theta*(L+mu)/(L+1)/(mu+1))^2;
else
    WC   = (2-theta)/4/mu/(L^2+1) * ...
        (theta*(1-2*mu+L^2)-2*mu*(L^2-1))*...
        (theta*(1+2*mu+L^2)-2*(mu+1)*(L^2+1))/...
        (theta*(1+2*mu-L^2)-2*(mu+1)*(1-L^2));
end
WC
