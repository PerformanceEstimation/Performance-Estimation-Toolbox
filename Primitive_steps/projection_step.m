function x = projection_step(x0,ind)

g_ind=Point('Point');
x=x0-g_ind;
f=Point('Function value');
ind.AddComponent(x,g_ind,f);

end

