function q=qf(p,u)
Kx=kron([[1,0,0];[0,1,0];[0,0,1]],p.mat.Kx); par=u(p.nu+1:end); uold=p.u(1:p.nu);
uox=Kx*uold;q=uox'*(u(1:p.nu)-uold); 
end
