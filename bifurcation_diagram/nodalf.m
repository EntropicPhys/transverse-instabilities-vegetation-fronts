function f=nodalf(p,u)
par=u(1:p.nu+1:end); b=u(1:p.np); w=u(p.np+1:2*p.np); h=u(2*p.np+1:3*p.np);
pp=par(1); Lam=par(2); nw=par(3); nh=par(4); ee=par(5); Rw=par(6); Rh=par(7);
Ga=par(8); A=par(9); Q=par(10); ff=par(11); k=par(12); m=par(13);
Lw=nw./(1+Rw*b/k); Lh=nh./(1+Rh*b/k); I=A*(b+Q*ff)./(b+Q);

f1=Lam*w.*b(1-b/k).*(1+ee*b).^2-m*b;
f2=I.*h - Lw.*w - Ga*w.*b.*(1+ee*b).^2;
f3=pp-I.*h-Lh.*h;
f=[f1;f2;f3];
