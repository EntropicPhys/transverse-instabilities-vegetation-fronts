function p=bwinit(p,lx,nx,par,b0,w0,h0,dir)
p=stanparam(p); %dir='init';
if ~exist(['./' dir],'dir'); mkdir(['./' dir]); end
p=setfn(p,dir); screenlayout(p);
p.nc.neq=3; p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs_ref;
p.fuha.outfu=@sgbra; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; p.vol=2*lx;
p.np=pde.grid.nPoints; p.nu=p.nc.neq*p.np; p.sol.xi=1/(p.nu);
p=setfemops(p); p.nc.ilam=1; p.nc.lammin=-10; p.nc.lammax=1000; p.nc.nsteps=20; 
b=b0*ones(p.np,1); w=w0*ones(p.np,1); h=h0*ones(p.np,1); p.u=[b;w;h;par];
p.sol.ds=0.01; % p.sw.verb=2;
p.plot.pbcmp=2; p.plot.cl={'black','blue','red'}; plotsol(p);
