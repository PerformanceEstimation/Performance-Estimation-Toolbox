function A_Minimizer_of_a_sum
% This example shows how the PESTO toolbox can also be used to explore
% properties of (strongly) convex functions that are not directly related
% to optimization algorithms:
% we consider the following problem: Given two functions f, g both
% mu-strongly convexe and L-smooth, how far can the minimizer xsfg of (f+g) 
% be from the mean of the minimizers xsf and xsg of respectively f and g. 
%
% One can verify that the answer to this problem 
% (i) is invariant if mu and L are multiplied by a same constant, and hence
% only depends on the condition number kappa = L/mu
% (ii) scales linearly with ||xsf-xsg||, which we will thus assume to be 1
%
% The code computes the worst-case distance for various values of kappa and
% plots its dependence on kappa.


% number of tests and vector of kappa to be tested
n_test = 80;
kappa = 1.025.^(1:n_test);
verbose = 0; % verbose mode of the toolbox


for k = 1:n_test;
    
    % (0) Information about current test    
    disp('---------------------------------------------------')
    disp(['Test ', num2str(k) , ' of ', num2str(n_test)])
    disp(['condition number = ' , num2str(round(kappa(k)*10^3)/10^3)])
    disp('---------------------------------------------------')

    % (1) Declaration of the "PEP"
    P=pep();

    % (2) functions parameters (using scale_invariance)
    param.mu=1;	% Strong convexity parameter
    param.L=kappa(k);      % Smoothness parameter

    % (3) functions declarations
    f= P.DeclareFunction('SmoothStronglyConvex',param);
    g= P.DeclareFunction('SmoothStronglyConvex',param);
    fg = f+g;
    
    % (4) decleration of the minimizers
    [xsf,fs] = f.OptimalPoint();
    [xsg,gs] = g.OptimalPoint();
    [xsfg,fgs] = fg.OptimalPoint();
     
    % (5) Contraint on the minimizer (have to be declared as
        P.AddConstraint((xsf-xsg)^2<=1);
    % note: the problem will naturally force (xsf-xsg)^2 = 1
    
    % (6) Criterion to be maximized
    P.PerformanceMetric((xsfg-(xsf+xsg)/2 )^2 ); 
    
    % (7) Solving the Pep
    P.solve(verbose);

    % (8) Evaluation and storage the output
    distance_sq(k) =double((xsfg-(xsf+xsg)/2 )^2); 
 
end


% representation of the results
distance = sqrt(distance_sq)

figure
plot(kappa,distance_sq)
xlabel('condition number');
ylabel('||x^s_{fg} - \frac{1}{2}(x^s_f+x^s_g)||^2')
title('square distance to average minimizers f, g')

figure
plot(kappa,distance)
xlabel('condition number');
ylabel('||x^s_{fg} - \frac{1}{2}(x^s_f+x^s_g)||')
title('distance to average minimizers f, g')


disp('************************************************************')
disp('Computation ended')
disp('The result should show the square distance is asymptotically')
disp('linear in the condition number')
disp('************************************************************')
end