function cons=Interpolation_ConvexSupport(pt1,pt2,D,R)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=(pt1.x*(pt2.g-pt1.g)<=0);
    if D~=Inf
        cons=cons+((pt1.g-pt2.g)^2-D^2<=0);
    end
else
    cons=(pt1.x*pt1.g-pt1.f==0);
    if R~=Inf
        cons=cons+((pt1.g)^2-R^2<=0);
    end
end
end