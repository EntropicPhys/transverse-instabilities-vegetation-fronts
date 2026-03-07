function qu=qfder(p,u)
Kx=kron([[1,0,0];[0,1,0];[0,0,1]],p.mat.Kx);qu=(Kx*p.u(1:p.nu))';
end
