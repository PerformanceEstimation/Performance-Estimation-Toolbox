function cons=Interpolation_ConvexSupport(pt1,pt2,M1,M2)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=(pt1.x*(pt2.g-pt1.g)<=0);
    if M1~=Inf
        cons=cons+((pt1.g-pt2.g)^2-M1^2<=0);
    end
else
    cons=(pt1.x*pt1.g-pt1.f==0);
    if M2~=Inf
        cons=cons+((pt1.g)^2-M2^2<=0);
    end
end
end