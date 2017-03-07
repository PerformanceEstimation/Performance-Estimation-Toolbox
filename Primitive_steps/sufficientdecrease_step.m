function [ x ] = sufficientdecrease_step(x0,func,expr)
x=Point('Point');
[g,f]=func.oracle(x);
[g0,f0]=func.oracle(x0);
func.AddConstraint(f-f0+expr);

end

