function [p,t1,vel]=tintfreeze(p,t0,dt,nt,pmod,vel,vmod)
n=0;t=t0;par=p.u(p.nu+1:end);dw=par(11);dh=par(12);K=p.mat.K;Kx=p.mat.Kx;M=p.mat.M;
bK=kron([[1,0,0];[0,dw,0];[0,0,dh]],K);
Lam=M+dt*bK; [L,U,P,Q,R]=lu(Lam);
while n<nt
      f=nodalf(p,p.u);
      G=bK*p.u(1:p.nu)-M*f;
      KX=kron([[1,0,0];[0,1,0];[0,0,1]],Kx);
      ux=KX*p.u(1:p.nu);cs=(ux')*G/((ux')*ux);
      g=M*p.u(1:p.nu)+dt*(M*f+cs*KX*p.u(1:p.nu));
      if (mod(n,vmod)==0) 
          vel=[vel [t;cs]];
      end
      p.u(1:p.nu)=Q*(U\(L\(P*(R\g))));
      n=n+1;
      t=t+dt;
      if (mod(n,pmod)==0)
          plotsol(p,p.plot.ifig, p.plot.pcmp, p.plot.pstyle)
      end
end
t1=t;
