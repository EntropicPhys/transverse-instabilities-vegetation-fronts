%% bwh co-moving front cont. 
format compact; keep pphome; 
%% dim. params 
PP=144.5;Lam=0.33/2;Nw=6.4;Nh=7;EE=1.55;Rw=0.73;Rh=0.3;
Ga=1.9;A=20;Q=0.065;ff=0.4;K=0.8;M=3.55;Db=1;Dw=90;Dh=3000;
%% non-dim params. 
nw=Nw/M;nh=Nh/M;a=A/M;q=Q/K;et=EE*Kl;
pp=PP*Lam/M^2;ga=Ga*K/M;dw=Dw/Db;dh=Dh/Db;s=0;
par=[pp;nw;nh;et;Rw;Rh;ga;a;q;ff;dw;dh;s];
%% print init pos. phase diag.
names={'pp','nw','nh','et','Rw','Rh','ga','a','q','ff','dw','dh','s'};
for i=1:length(par);fprintf('%s=%.6g\n',names{i},par(i));end
%% frozen state.
dir= ['cont/dw' num2str(i) '/eta' num2str(j)];
lx=360;nx=1000; p=[]; p=bwhinit(p,lx,nx,par,0.455,0.845,1,dir);
nt=1e4;pmod=350;vmod=100;dt=1e-2;vel=[];t0=0;
[p,t1,vel]=tintfreeze(p,t0,dt,nt,pmod,vel,vmod);
%% cont.
p.nc.ilam=[4 13]; p.nc.nq=1; p.u(p.nu+13)=vel(2,end);p.sw.bifcheck=1;
p.plot.bpcmp=13;p.nc.dsmin=1e-5; p.plot.pstyle=5; p.fuha.qf=@qf;
p.fuha.qfder=@qfder; p.sw.qjac=1; p.sw.spcalc=0; p.nc.dsmax=3e-2;
p.sol.ds=p.nc.dsmax;p=cont(p,nc);
