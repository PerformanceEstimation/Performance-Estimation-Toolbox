function cons=Interpolation_Convex(pt1,pt2)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=((pt1.f-pt2.f+pt1.g*(pt2.x-pt1.x))<=0);
else
    cons=[];
end

end