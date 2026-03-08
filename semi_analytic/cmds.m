%% bwh co-moving front cont. 
format compact; keep pphome; 
%% dim. params 
PP=144.5;Lam=0.33/2;Nw=6.4;Nh=7;EE=1.55;Rw=0.73;Rh=0.3;
Ga=1.9;A=20;Q=0.065;ff=0.7;K=0.8;M=3.55;Db=1;Dw=90;Dh=3000;
%% non-dim params. 
nw=Nw/M;nh=Nh/M;a=A/M;q=Q/K;et=EE*K;
pp=PP*Lam/M^2;ga=Ga*K/M;dw=Dw/Db;dh=Dh/Db;s=0;
par=[pp;nw;nh;et;Rw;Rh;ga;a;q;ff;dw;dh;s];
%% frozen state.
dir= strrep(num2str(ff),'.','p');
lx=360;nx=1500; p=[]; 
p=bwhinit(p,lx,nx,par,dir);
%%
nt=1e3;pmod=350;vmod=100;dt=1e-2;vel=[];t0=0;
[p,t1,vel]=tintfreeze(p,t0,dt,nt,pmod,vel,vmod);
%% cont.
p.nc.ilam=[4 13]; p.nc.nq=1; p.u(p.nu+13)=vel(2,end);p.sw.bifcheck=1;
p.plot.bpcmp=13;p.nc.dsmin=1e-5; p.plot.pstyle=5; p.fuha.qf=@qf;
p.fuha.qfder=@qfder; p.sw.qjac=1; p.sw.spcalc=0; p.nc.dsmax=3e-2;
p.file.smod=1;p.sol.ds=p.nc.dsmax;p=cont(p,20);
%% trans. stability.
indx = sort(getlabs(dir));lambda=[];eta=[];eig_trans=[];eig_adj=[];
for k = 1:length(indx)
    pt = ['pt' mat2str(indx(k))];
    p0 = loadp(dir, pt);
    [lambda0,eta0,eig_trans0,eig_adj0]=compute_lam_eta(p0);
    lambda = [lambda; lambda0 ];eta=[eta;eta0];
    eig_trans=[eig_trans; eig_trans0];
    eig_adj=[eig_adj;eig_adj0];
end
%% plot stabilty curve
figure(11)
plot(eta,lambda,'.-');
title('spectrum slope');grid on
xlabel('$\eta$','Interpreter','latex');ylabel('$\lambda_1$','Interpreter','latex');
figure(12)
subplot(2,1,1)
plot(eta,eig_trans,'.-')
title('eigenvalue translational mode');grid on
xlabel('$\eta$','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');
subplot(2,1,2)
plot(eta,eig_adj,'.-')
title('Eigenvalue kernel adjoint');grid on
xlabel('$\eta$','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');
