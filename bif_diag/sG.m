function r=sG(p,u)
par=u(p.nu+1:end); dw=par(14);dh=par(15); 
K=p.mat.K; M=p.mat.M; bK=kron([[1,0,0];[0,dw,0];[0,0,dh]],K);
r=bK*u(1:p.nu) - M*nodalf(p,u); 
