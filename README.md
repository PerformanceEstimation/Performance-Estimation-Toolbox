# Performance-Estimation-Toolbox (PESTO)

This code comes jointly with the following [`reference`](PESTO_CDC2017_FINAL.pdf).

> [1] Taylor, Adrien B., Julien M. Hendrickx, and François Glineur. "Performance Estimation Toolbox (PESTO): automated worst-case analysis of first-order optimization methods." Proceedings of the 56th IEEE Conference on Decision and Control (CDC 2017).

Date:    Nov 2019

Version: Nov 2019

#### Authors
- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)
- [**Julien Hendrickx**](https://perso.uclouvain.be/julien.hendrickx/index.html)
- [**François Glineur**](https://perso.uclouvain.be/francois.glineur/)

#### Acknowledgments

The authors would like to thank [**Francois Gonze**](https://perso.uclouvain.be/francois.gonze/) from UCLouvain and **Yoel Drori** from Google Inc. for their feedbacks on preliminary versions of the toolbox.

Additional material was incorporated thanks to:
- [**Ernest Ryu**](https://www.math.ucla.edu/~eryu/) (UCLA), [**Carolina Bergeling**](http://www.control.lth.se/personnel/carolina-bergeling.html) (Lund), and [**Pontus Giselsson**](http://www.control.lth.se/personnel/pontus-giselsson/) (Lund) (monotone operators and splitting methods),
- [**Francis Bach**](https://www.di.ens.fr/~fbach/) (Inria & ENS Paris) (stochastic methods and potential functions),
- **Radu-Alexandru Dragomir** (ENS Paris & TSE), [**Alexandre d’Aspremont**](https://www.di.ens.fr/~aspremon/) (CNRS & ENS Paris), and [**Jérôme Bolte**](https://www.tse-fr.eu/people/jerome-bolte) (TSE) (Bregman divergences and new notions of smoothness).

Last but not least, we thank **Loic Estève** (Inria) for technical support.

#### Performance Estimation Framework

The toolbox implements the performance estimation approach as developped in the following articles:
 > [2] Taylor, Adrien B., Julien M. Hendrickx, and François Glineur. "Smooth strongly convex interpolation and exact worst-case performance of first-order methods." Mathematical Programming 161.1-2 (2017): 307-345.

 > [3] Taylor, Adrien B., Julien M. Hendrickx, and François Glineur. "Exact worst-case performance of first-order methods for composite convex optimization." SIAM Journal on Optimization 27.3 (2017): 1283-1313

Note that the approach of using semidefinite programming for obtaining worst-case guarantees was originally introduced in

 > [4] Drori, Yoel, and Marc Teboulle. "Performance of first-order methods for smooth convex minimization: a novel approach." Mathematical Programming 145.1-2 (2014): 451-482

## Introduction to the toolbox

The document [`PESTO_CDC2017_FINAL`](PESTO_CDC2017_FINAL.pdf) contains the reference paper for the toolbox.
This paper contains a simplified and quick general introduction to the theory underlying the toolbox, and to its use.

The document [`UserGuide`](UserGuide.pdf) contains more detailed information and examples about the use of the toolbox.

The general purpose of the toolbox is to help the researchers producing worst-case guarantees for their favorite first-order methods.

## Setup

**Note:** This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).

Once YALMIP and the SDP solver installed (type 'yalmiptest' for checking the installation went fine); the toolbox can simply be installed by typing
```Matlab
Install_PESTO
```
in the Matlab prompt (which only adds the required folders to your Matlab path).
You can now execute the demo files for a step by step introduction to the toolbox.


## Example

The folder [Examples](/Examples) contains numerous introductory examples to the toolbox.

The files [`demo`](demo.m) summarizes the currently available step-by-step demonstration files and their targets. As an example, the first demonstration file can be executed by typing
```Matlab
demo1
```
in the prompt.


Among the other examples, the following code (see [`GradientMethod`](/Examples/GradientMethod.m)) generates a worst-case scenario for the gradient method applied to a smooth (possibly strongly) convex function.

```Matlab
% In this example, we use a fixed-step gradient method for
% solving the L-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x).
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the gradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.mu=0;	% Strong convexity parameter
param.L=1;      % Smoothness parameter

F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		  % x0 is some starting point
[xs,fs]=F.OptimalPoint(); 	  % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1); % Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm
h=1/param.L;		% step size
N=10;		% number of iterations

x=x0;
for i=1:N
    x=x-h*F.gradient(x);
end
xN=x;

% (4) Set up the performance measure
fN=F.value(xN);                % g=grad F(x), f=F(x)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(fN-fs)   % worst-case objective function accuracy

% The result should be
% param.L/2/(2*N+1)

```



 


