function cons=Interpolation_ConvexIndicator(pt1,pt2,M1,M2)

if ~(pt1.x.isEqual(pt2.x))
    cons=(pt1.g*(pt2.x-pt1.x)<=0);
    if M1~=Inf
        cons=cons+((pt1.x-pt2.x)^2-M1^2<=0);
    end
else
    cons=(pt1.f==0);
    if M2~=Inf
        cons=cons+((pt1.x)^2-M2^2<=0);
    end
end

end