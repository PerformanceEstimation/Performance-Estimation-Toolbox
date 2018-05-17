function plotSmoothInterp(F,L)
% TODO: faire avec mu aussi
% TODO changer longueur des tangentes
if nargin<2
    L=1;
end
tags = F.GetTags();
n = length(tags);
if n == 0
    fprintf('No tags found ; only tagged points can be plotted.\n');
    return;
end
dim = length(double(F.gradient(tags{1})));
if dim ~= 1
   fprintf('Plot requires a 1D function ; points here have dimension %d\n', dim);
   return;
end
for i=1:n
    [gobj,fobj,xobj] = F.oracle(tags{i});
    x(i) = double(xobj);
    g(i) = double(gobj);
    f(i) = double(fobj);
end
[x,ord] = sort(x);
f = f(ord);
g = g(ord);
tags = tags(ord);

title(sprintf('Smooth convex interpolated function (L=%d) with %d points', L, n));
hold on;
xmin = min(x); xmax = max(x); xrange = xmax-xmin;
fmin = min(f); fmax = max(f); frange = fmax-fmin;
zoomfact = 0.1;
axis([xmin-xrange*zoomfact xmax+xrange*zoomfact fmin-frange*zoomfact fmax+frange*zoomfact]);
tangentdelta = xrange * 0.02;
textfact = 0.025;
for i=1:n
    plot(x(i),f(i),'r.','MarkerSize',25);
    text(x(i)+xrange*textfact,f(i)+frange*textfact,tags(i));
    plot(x(i)+tangentdelta*[-1 1], f(i)+tangentdelta*g(i)*[-1 1],'b:','LineWidth',2); 
end
for i=1:n-1
    len2 = (g(i+1)-g(i))/L;
    ginc = L*len2^2/2;
    otherinc = f(i+1)-f(i)-g(i)*(x(i+1)-x(i)) - ginc;
    if abs(otherinc) < 1e-8
        len3 = 0;
    else
        len3 = (f(i+1)-f(i)-g(i)*(x(i+1)-x(i))-ginc)/g(i+1);
    end
    len1 = x(i+1)-x(i)-len3-len2;
    if len1>1e-4
        plot([x(i) x(i)+len1],[f(i) f(i)+len1*g(i)],'k-','LineWidth',2);    
        plot(x(i)+len1, f(i)+len1*g(i),'x','MarkerSize',15);
    end
    if len3>1e-4
        plot([x(i+1)-len3 x(i+1)],[f(i+1)-len3*g(i+1) f(i+1)],'k-','LineWidth',2);
        plot(x(i+1)-len3, f(i+1)-len3*g(i+1),'x','MarkerSize',15);
    end
    if len2>1e-4
        mid = linspace(x(i)+len1,x(i+1)-len3);
        fmid = f(i)+len1*g(i) + (mid-x(i)-len1)*g(i) + (mid-x(i)-len1).^2/2/L;
        plot(mid,fmid,'k-','LineWidth',2);
    end
    a=1;    
end
    
    
    %[x(i+1)-x(i) f(i+1)-f(i) (g(i+1)+g(i))/2*(x(i+1)-x(i))]
hold off;
figure(gcf);
return;    
