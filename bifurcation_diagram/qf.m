function q=qf(p,u)
Kx=p.mat.Kx; uold=p.u(1:p.nu); uox=Kx*uold; q=uox'*(u(1:p.nu)-uold);
end
