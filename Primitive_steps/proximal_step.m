function [x] = proximal_step(x0,func,gamma)

g_imp=Point('Point');
x=x0-gamma*g_imp;
f=Point('Function value');
func.AddComponent(x,g_imp,f);

end

