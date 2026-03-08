%% bwh model 
close all; format compact; keep pphome; 
%% dim. params
PP=200;Lam=0.33/2;Nw=6.4;Nh=7;EE=1.355;Rw=0.73;Rh=0.3;Ga=1.9;
A=20;Q=0.065;ff=0.4;K=0.8;M=3.55;Db=1;Dw=90;Dh=3000; 
par=[PP; Lam; Nw; Nh; EE; Rw; Rh; Ga; A; Q; ff; K; M; Dw; Dh];
%% dir. & init. cond.
lx=100;nx=600;b0=0.622;w0=28.5;h0=PP/(A*(b0+Q*ff)/(b0+Q)+Nh/(1+Rh*b0/K));
p=[]; p=bwhinit(p,lx,nx,par,b0,w0,h0,'fig2/hom');
%% hom. state cont. 
p.nc.tol=1e-8;p.nc.dsmax=2;p.sol.ds=-1e-1;p.nc.neig=3;p=cont(p,140);
%% first pp cont.
aux=[];aux.m=3; aux.besw=0; p0=cswibra('fig2/hom','bpt1',aux);
p0.nc.tol=1e-10;p0.nc.neig=5;p0.sw.bifcheck=2;
p0.nc.dsmin=1e-6; p0.nc.dsmax=2; dirT='fig2/T1';p0.nc.eigref=-1;
p=gentau(p0,[1]);p=setfn(p,dirT); p.sol.ds=1e-3; p=cont(p,160);
%% bare-soil cont. 
PP0=230;p.u(p.nu+1)=PP0; b0=0; w0=A*ff*PP0/(Nw*A*ff+Nw*Nh); h0=PP0/(A*ff+Nh);
p=resetc(p);p.u(1:p.np)=zeros(p.np,1);p.u(p.np+1:2*p.np)=w0*ones(p.np,1);p.u(2*p.np+1:3*p.np)=h0*ones(p.np,1);
p=setfn(p,'fig2/hom1');p.nc.tol=1e-8;p.nc.dsmax=10;p.sol.ds=6e-1;p=cont(p,60);
%% plot sol. branches
plotbra('fig2/hom','cl','#408000')
plotbra('fig2/T1','cl','r')
plotbra('fig2/hom1','cl','#FF8000')
axis([110 300 0 0.65]); xticks([110 200 280]); yticks([0 0.2 0.4 0.6])
ylabel('$||B|| [kg/m^2]$','Interpreter','latex')
xlabel('$P [mm/y]$','Interpreter','latex')
box on 
