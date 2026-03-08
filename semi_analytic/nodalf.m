function f=nodalf(p,u)
par=u(p.nu+1:end); n=p.np; b=u(1:n); w=u(n+1:2*n); h=u(2*n+1:3*n);
pp=par(1); nw=par(2); nh=par(3); et=par(4); Rw=par(5); Rh=par(6);
ga=par(7); a=par(8); q=par(9); ff=par(10); 
I=a*(b+q*ff)./(b+q); Lw=nw./(1+Rw*b); Lh=nh./(1+Rh*b); Gw=ga*b.*(1+et*b).^2;
f1=w.*b.*(1-b).*(1+et*b).^2-b;
f2=I.*h-Lw.*w-Gw.*w;
f3=pp-I.*h-Lh.*h;
f=[f1;f2;f3];
end
