

Nmax=30; eps_values=[.1 .2 .5];
N=1:Nmax;
performance=zeros(length(eps_values)+1,Nmax);
performance(1,:)=2./(N.^2+5*N+6);%eps=0:
for i=N
    for j=1:length(eps_values)
        fprintf('Case N:%d/%d, eps:%d/%d\n',i,Nmax,j,length(eps_values));
        performance(j+1,i)=CDC_fastgradientmethod(i,eps_values(j));
    end
end

%%

loglog(N,performance(1,:),'-g','linewidth',2); set(gca,'FontSize',14); hold on;
loglog(N,performance(2,:),'-r','linewidth',2); 
loglog(N,performance(3,:),'-b','linewidth',2);
loglog(N(1:15),performance(4,1:15),'-m','linewidth',2);
xlabel('Iteration count (log scale)','Fontsize',14);
ylabel('Objective function accuracy (log scale)','Fontsize',14);
% legend('N=2','N=5','N=10');
print -depsc InexactFGM_WC.eps
