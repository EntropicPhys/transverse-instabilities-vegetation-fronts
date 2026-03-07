function r=sG(p,u)
nu=p.nu;par(nu+1:end);dw=par(11);dh=par(12);s=par(13);
K=p.mat.K;Kx=p.mat.Kx;M=p.mat.M;
bK=kron([[1,0,0];[0,dw,0];[0,0,dh]],K);sKx=kron([[s,0,0];[0,s,0];[0,0,s]],Kx);
r=(bK-sKx)*u(1:nu)-M*nodalf(p,u);
