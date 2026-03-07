function r=sG(p,u)
par=u(p.nu+1:end); dw=par(14); dh=par(15); s=par(16); 
K=p.mat.K; Kx=p.mat.Kx; bK=kron([1,0,0];[0,dw,0];[0,0,dh],K);
M=p.mat.M; r=(bK-s*Kx)*u(1:p.nu) - M*nodalf(p,u);
