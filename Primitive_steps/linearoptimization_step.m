function [x] = linearoptimization_step(dir,ind)

g_imp=-dir;
x=Point('Point');
feas=Point('Function value');
ind.AddComponent(x,g_imp,feas);

end

