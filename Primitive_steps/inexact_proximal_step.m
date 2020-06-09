function [x,gx,fx,w,v,fw,epsVar] = inexact_proximal_step(x0,func,gamma,opt)

% [x] = inexact_proximal_step(x0,func,gamma,parameters)
%
% This routine performs an inexact proximal step with step size gamma,
% starting from x0, and on function func. That is, it performs:
%       x = x0-gamma*(v+e), where g is a (epsilon-sub)gradient of func at
%                           x, and e is some computation error whose 
%                           characteristics are provided in the "settings"
%                           structure.
%
% Input: - starting point x0
%        - function func on which the (sub)gradient will be evaluated
%        - step size gamma of the proximal step
%        - inaccuracy parameters in the "opt" structure.
%
% Output: - x = x0-gamma*(v-e), where v is a eps-subgradient of func at x.
%         - gx a subgradient of func at x
%         - fx function func evaluated at x
%         - w a primal point (possibly = x) such that v is a subgradient of
%         func at w.
%         - v is a subgradient of func at w.
%         - fw function func evaluated at w 
%         - epsVar is the required accuracy (a variable), which the user
%         can (should) bound.
%   NOTE: v is an epsilon-subgradient of func at x. In order to compute
%   this epsilon, we simply have to do: epsilon = fx - fw - v*(x-w).
%
% Details on the 'opt' structure. It might contain:
%   - opt.criterion: see below.
%
% Details for opt.criterion:
%   - 'PD_gapI':
%           PD gap(x,v;x0) <= epsVar for the proximal subproblem.
%   - 'PD_gapII':
%           PD gap(x,v;x0) <= epsVar for the proximal subproblem,
%           with v \in \partial func(x).
%   - 'PD_gapIII':
%           PD gap(x,v;x0) <= epsVar for the proximal subproblem,
%           with v = (x_0-x)/gamma.
%   - 'Orip-style' (see [1] below)
%        Approximate proximal operator outputs x such that 
%        <v; e> + epsilon/gamma <= epsVar for the proximal subproblem,
%        with x = x_0 - gamma * ( v - e ) with v an epsilon-subgradient of
%        func at x.
%
%
%
% ORIP: optimized relatively inexact proximal point algorithm (see [1])
% [1] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators

switch opt.criterion
    case 'PD_gapI'
        % Approximate proximal operator outputs x such that 
        % PD gap(x,v;x0) <= epsVar
        % with v some dual variable, and
        % PD gap(x,v;x0)= 1/2*(x-x0+gamma*v)^2 + gamma *(func(x)+func*(v)-<v;x>),
        % in which we conveniently use fync*(v) = <v,w>-func(w) for some w such
        % that v\in\partial func(w).
                
        v   = Point('Point');
        w   = Point('Point');
        fw  = Point('Function value');
        func.AddComponent(w,v,fw);
        
        gx  = Point('Point');
        x   = Point('Point');
        fx  = Point('Function value');
        func.AddComponent(x,gx,fx);
        
        
        epsVar  = Point('Function value');
        
        e       = x-x0+gamma*v;
        eps_sub = fx - fw - v*(x-w);
        gap     = 1/2*e^2 + gamma*eps_sub;
        
        func.AddConstraint(gap<=epsVar);     
    case 'PD_gapII'
        % Approximate proximal operator outputs x such that 
        % ||e|| <= epsVar
        % with x = x_0 - gamma * ( gx - e ) and gx a subgradient of func at x.
                
        e   = Point('Point');
        gx  = Point('Point'); v = gx;
        x   = x0 - gamma * ( v - e ); w = x;
        fx  = Point('Function value'); fw = fx;
        func.AddComponent(x,gx,fx);
        epsVar  = Point('Function value');
                
        func.AddConstraint(e^2<= epsVar); 
    case 'PD_gapIII'
        % Approximate proximal operator outputs x such that 
        % gamma * (fx - fw - v*(x-w))  <= epsVar
        % v=(x0-x)/gamma and w a point such that v a subgradient of func at w.
                
        x   = Point('Point');
        v   = (x0-x)/gamma;
        w   = Point('Point');
        fw  = Point('Function value');
        gx  = Point('Point');
        fx  = Point('Function value');
        func.AddComponent(w,v,fw);
        func.AddComponent(x,gx,fx);
        epsVar  = Point('Function value');
                
        eps_sub = fx - fw - v*(x-w);
        
        func.AddConstraint(gamma*eps_sub<=epsVar);       
        
    case 'Orip-style'
        % Approximate proximal operator outputs x such that 
        % <v; e> + epsilon/gamma <= epsVar
        % with x = x_0 - gamma * ( v - e ) with v an epsilon-subgradient of
        % func at x.
        
        v   = Point('Point');
        w   = Point('Point');
        fw  = Point('Function value');
        func.AddComponent(w,v,fw);
        
        e   = Point('Point');
        gx  = Point('Point');
        x   = x0 - gamma * ( v - e );
        fx  = Point('Function value');
        func.AddComponent(x,gx,fx);
        epsVar  = Point('Function value');
                
        eps_sub = fx - fw - v*(x-w);
        gap     = e*v + eps_sub/gamma;
        
        func.AddConstraint(gap<=epsVar);        
    otherwise
        fprintf('No criterion chosen for inexact prox\n');
end
end 

