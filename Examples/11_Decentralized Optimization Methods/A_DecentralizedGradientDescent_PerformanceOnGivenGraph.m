clear all; clc;
% In this example, we consider K iterations of the decentralized subgradient
% descent with N agents that each holds a local convex function Fi with bounded subgradients
% for solving the following decentralized problem:
%   min_x F(x);     where F(x) is the sum of local functions Fi.
% Agents communicate through a given communication network, represented by the communication matrix W.
%
% This script calls the function DGD_exact_perf to compute the exact worst-case performance of DGD,
% with repsect to the performance measure  F(xav)-F(xs) where xav the average of all the iterates (for each iteration and each agent).
% The initial iterates satisfy ||x0 - x*||^2 <= IC^2, for all agents.
%
% For details, see
%   Colla, Sebastien, and Julien M. Hendrickx. "Automated Worst-Case
%   Performance Analysis of Decentralized Gradient Descent." (2021)

K = 10;                 % Number of iterations of DGD
alpha = 1./sqrt(K);     % Step-size used in DGD (constant)
%alpha = 1./(1:K);       % Alternative: Step-sizes used in DGD (diminishing)
N = 3;                  % Number of agents
W = 1/N*ones(N,N);      % Communication matrix
IC = 1;                 % Constant for the initial condition: ||x0 - xs||^2 <= IC^2
equalStart = 1;         % All agents starts with the same iterate x0
fctClass = 'ConvexBoundedGradient'; % Class of functions to consider for the worst-case
fctParam.R = 1;         % Bounded subgradient constant ||g||^2 <= R^2.
avgAll = 1;             % The performance bound considers the average iterates 'xav' over all agents and all iterates: F(xav) - F(xs).
verbose = 1;            % Print the problem and the results

[wc, out] = DGD_exact_perf(K,alpha,N,W,IC,equalStart,fctClass,fctParam,avgAll,verbose);