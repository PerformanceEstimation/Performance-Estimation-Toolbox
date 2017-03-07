function cons=Interpolation_StronglyConvexBoundedDomain(pt1,pt2,M1,M2,mu)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+mu/2*(pt1.x-pt2.x)^2)<=0);
    if M1~=Inf && pt1.num<pt2.num
        cons=cons+((pt1.x-pt2.x)^2-M1^2<=0);
    end
else
    if M2~=Inf
        cons=((pt1.x)^2-M2^2<=0);
    else
        cons=[];
    end
end

end