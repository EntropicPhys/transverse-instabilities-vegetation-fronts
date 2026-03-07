function [lam,eta]=compute_lam_eta(p0)
Mblk=p0.mat.M;par=p0.u(p0.nu+1:end);dw=par(11);dh=par(12);eta=par(4);n=p0.np;xx=getpte(p0);
[Vt,Dt]=eigs(LOp(p0,0),Mblk,3,-0.1);[Vad,Dad]=eigs(LAp(p0),Mblk,3,-0.1);
F1=Vad(1:n,1).*Vt(1:n,1)+dw*Vad(n+1:2*n,1).*Vt(n+1:2*n,1)+dh*Vad(2*n+1:3*n,1).*Vt(2*n+1:3*n,1);
F2=Vad(1:n,1).*Vt(1:n,1)+Vad(n+1:2*n,1).*Vt(n+1:2*n,1)+Vad(2*n+1:3*n,1).*Vt(2*n+1:3*n);
lam=-trapz(xx,F1)/trapz(xx,F2);
end
