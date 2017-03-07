function [x] = gradient_step(x0,func,gamma)

[g,~]=func.oracle(x0);
x=x0-gamma*g;

end

