function p=bwhinit(p,lx,nx,par,b0,w0,h0,dir) % Initialization of continuation and FEM operators/space based on pdetoolbox  
p=stanparam(p); % sets main continuation parameter and functions 
p=setfn(p,dir);  % Creates dir as directory.
p.fuha.outfu=@sgbra; % sgbra.m is set as file to compute continuation outputs.

% PDE objects and key quantities are defined; see stanparam for details.
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; p.vol=2*lx;
p.np=pde.grid.nPoints; p.nu=p.nc.neq*p.np; p.sol.xi=1/(p.nu);
p=setfemops(p); p.nc.ilam=1; p.nc.lammin=0;  
p.nc.neq=3; p.sw.sfem=-1; 

b=b0*ones(p.np,1); w=w0*ones(p.np,1); h=h0*ones(p.np,1); p.u=[b;w;h;par]; % initial conditions and assembly of p.u
p.sol.ds=0.01; p.sw.verb=2; % initial arc-length step size and terminal print protocol
screenlayout(p); p.plot.bpcmp=2; p.plot.cl={'black','blue','red'}; plotsol(p); % generic plot parameters; see staparam for details.
