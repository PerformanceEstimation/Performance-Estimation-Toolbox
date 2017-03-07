function cons=Interpolation_SmoothConvexBoundedGradient(pt1,pt2,M1,M2,L)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    if L~=Inf
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x)+1/(2*L)*(pt1.g-pt2.g)^2)<=0);
    else
        cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x))<=0);
    end
    if M1~=Inf
        cons=cons+((pt1.g-pt2.g)^2-M1^2<=0);
    end
else
    if M2~=Inf
        cons=((pt1.g)^2-M2^2<=0);
    else
        cons=[];
    end
end

end