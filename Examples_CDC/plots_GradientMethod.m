
Nmax=2; nb_gamma=50;
gamma_vec=linspace(0,2,nb_gamma);
performance1=zeros(nb_gamma,1);
for i=1:nb_gamma
    fprintf('Case %d on %d (1/3)\n',i,nb_gamma);
    gamma=gamma_vec(i);
    performance1(i)=CDC_gradientmethod(Nmax,gamma*ones(Nmax,1));
end
Nmax=5; 
performance2=zeros(nb_gamma,1);
for i=1:nb_gamma
    fprintf('Case %d on %d (2/3)\n',i,nb_gamma);
    gamma=gamma_vec(i);
    performance2(i)=CDC_gradientmethod(Nmax,gamma*ones(Nmax,1));
end

Nmax=10; 
performance3=zeros(nb_gamma,1);
for i=1:nb_gamma
    fprintf('Case %d on %d (3/3)\n',i,nb_gamma);
    gamma=gamma_vec(i);
    performance3(i)=CDC_gradientmethod(Nmax,gamma*ones(Nmax,1));
end

plot(gamma_vec,performance1,'-g','linewidth',2); set(gca,'FontSize',14); hold on;
plot(gamma_vec,performance2,'-r','linewidth',2); 
plot(gamma_vec,performance3,'-b','linewidth',2);
xlabel('Step size','Fontsize',14);
ylabel('Residual gradient norm','Fontsize',14);
% legend('N=2','N=5','N=10');
print -depsc GradientWC.eps

