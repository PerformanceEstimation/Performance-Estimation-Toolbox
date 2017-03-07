function cons=Interpolation_Smooth(pt1,pt2,L)

if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    cons=((pt1.f-pt2.f-L/4*(pt1.x-pt2.x)^2+1/2*(pt1.g+pt2.g)*(pt1.x-pt2.x)+1/(4*L)*(pt1.g-pt2.g)^2)<=0);
else
    cons=[];
end

end