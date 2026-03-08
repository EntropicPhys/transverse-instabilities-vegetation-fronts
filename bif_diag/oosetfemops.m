function p=oosetfemops(p)
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1);
p.mat.K=K; p.mat.M=kron([[1,0,0];[0,1,0];[0,0,1]],M);
