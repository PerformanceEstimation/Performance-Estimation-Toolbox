function cons=Interpolation_SmoothStronglyConvex(pt1,pt2,mu,L)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    if L~=Inf
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+1/(2*(1-mu/L))*(1/L*(pt1.g-pt2.g)*(pt1.g-pt2.g).'+mu*(pt1.x-pt2.x)*(pt1.x-pt2.x).'-2*mu/L*(pt1.x-pt2.x)*(pt1.g-pt2.g).'))<=0);
    else
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+mu/2*(pt1.x-pt2.x)^2)<=0);
    end    
else
    cons=[];
end
end