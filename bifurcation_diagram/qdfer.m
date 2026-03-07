function qu=qdfer(p,u)
Kx=p.mat.Kx;qu=(Kx*p.u(1:p.nu))';
end
