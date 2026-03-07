function Gu=sGjac(p,u)
par=u(p.nu+1:end); pp=par(1); Lam=par(2); nw=par(3); nh=par(4); 
ee=par(5); Rw=par(6); Rh=par(7); Ga=par(8); A=par(9); Q=par(10); 
ff=par(11); k=par(12);m=par(13);dw=par(14);dh=par(15);s=par(16);
n=p.np;b=u(1:n);w=u(n+1:2*n);h=u(2*n+1:3*n); ov=ones(n,1);

f1b = -m*ov - Lam*w.*(1+ee*b).*(4*ee*b.^2-k+b*(2-3*ee*k))/k;
f1w = Lam*(1-b/k).*b.*(1+ee*b).^2;
f1h = 0*ov;
f2b = -A*(-1+ff)*Q*h./(b+Q).^2+k*nw*Rw*w./(k+b*Rw).^2 - Ga*(1+4*b*ee + ee^2*3*b.^2).*w;
f2w =  -k*nw./(k+b*Rw)-Ga*b.*(1+ee*b).^2 ;
f2h = A*(b+Q*ff)./(b+Q);
f3b = A*(-1+ff)*Q*h./(b+Q).^2 + k*nh*Rh*h./(k+Rh*b).^2;
f3w = 0*ov;
f3h = -A*(b+Q*ff)./(b+Q) - k*nh./(k+Rh*b);

Fu = [[spdiags(f1b,0,n,n), spdiags(f1w,0,n,n), spdiags(f1h,0,n,n)];
      [spdiags(f2b,0,n,n), spdiags(f2w,0,n,n), spdiags(f2h,0,n,n)];
      [spdiags(f3b,0,n,n), spdiags(f3w,0,n,n), spdiags(f3h,0,n,n)]];

K=p.mat.K; M=p.mat.M;Kx=p.mat.Kx;
bK=kron([[1,0,0];[0,dw,0];[0,0,dh]],K);%[[K 0*K 0*K];[0*K dw*K 0*K];[0*K 0*K dh*K]];
Gu=bK-s*Kx - M*Fu;
end
