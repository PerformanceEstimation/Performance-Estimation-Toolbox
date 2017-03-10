function cons=Interpolation_ConvexIndicator(pt1,pt2,D,R)

if ~(pt1.x.isEqual(pt2.x))
    cons=(pt1.g*(pt2.x-pt1.x)<=0);
    if D~=Inf
        cons=cons+((pt1.x-pt2.x)^2-D^2<=0);
    end
else
    cons=(pt1.f==0);
    if R~=Inf
        cons=cons+((pt1.x)^2-R^2<=0);
    end
end

end