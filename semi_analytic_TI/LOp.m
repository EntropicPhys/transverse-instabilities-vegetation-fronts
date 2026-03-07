function r=LOp(p,l)
par=p.u(p.nu+1:end);n=p.np;b=p.u(1:n);w=p.u(n+1:2*n);h=p.u(2*n+1:3*n);
pp=par(1); nw=par(2); nh=par(3); et=par(4); Rw=par(5); Rh=par(6); 
ga=par(7); a=par(8); q=par(9); ff=par(10); dw=par(11); dh=par(12); s=par(13);

I=a*(b+q*ff)./(b+q);dI=a*q*(1-ff)*ov./(b+q).^2;
Lw=nw*ov./(1+Rw*b);dLw=-nw*Rw*ov./(1+Rw*b).^2;
Lh=nh*ov./(1+Rh*b);dLh=-nh*Rh*ov./(1+Rh*b).^2;
Gw=ga*b.*(1+et*b).^2;dGw=ga*(1+et*b).*(1+3*et*b);

f1b=-1+w.*(1+3*et*(-2+et).*b.^2-4*(et^2)*b.^3+(4*et-2)*b);
f1w=(1-b).*b.*(1+et*b).^2;f1h=zeros(n,1);
f2b=dI.*h-dLw.*w-dGw.*w;f2w=-Gw-Lw;f2h=I;
f3b=-(dLh+dI).*h;f3w=zeros(n,1);f3h=-I-Lh;
ov=ones(n,1);
Lu = [[spdiags(ov,0,n,n), spdiags(0*ov,0,n,n), spdiags(0*ov,0,n,n)];
      [spdiags(0*ov,0,n,n), spdiags(dw*ov,0,n,n), spdiags(0*ov,0,n,n)];
      [spdiags(0*ov,0,n,n), spdiags(0*ov,0,n,n), spdiags(dh*ov,0,n,n)]];
K=p.mat.K;Kx=p.mat.Kx;M=p.mat.M;bK=kron([[1,0,0];[0,dw,0];[0,0,dh]],K);
sKx=kron([[s,0,0];[0,s,0];[0,0,s]],Kx);r=bK-sKx-M*Fu+l^2*M*Lu;
end
