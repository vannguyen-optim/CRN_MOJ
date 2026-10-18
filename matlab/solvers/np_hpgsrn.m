function out=np_hpgsrn(P,x,o)
% Reference implementation of HpgSRN, Algorithm 3.2, p=1/2.
% https://doi.org/10.1137/22M1482822
assert(strcmp(P.kind,'regression')&&isempty(P.C)&&all(isinf(P.lb)));
timer=tic;np_work_budget('start',o.max_qp);X=x;hist=[];elapsed=0;log={};
mu=1;xp=[];gp=[];status='iteration_budget';native=Inf;nlin=0;
for k=1:o.max_outer
 [~,R,ok]=np_certificate(P,x,[],o);if ok&&R<=o.tol,status='common_stop';break;end
 [fx,g,H]=np_eval(P,x);
 if ~isempty(xp)
  s=x-xp;r=g-gp;den=s'*s;
  if den>0,mu=min(1e20,max(1e-20,(s'*r)/den));end
 end
 xp=x;gp=g;accepted=false;
 for j=1:o.external_backtracks
  xb=np_half_prox(x-g/mu,P.alpha/mu);step=xb-x;
  if np_difference(P,x,xb)<=-o.hpg_alpha/2*(step'*step),accepted=true;break;end
  mu=10*mu;
 end
 if ~accepted,status='pg_backtracking_limit';break;end
 native=mu*norm(xb-x,Inf);
 % Algorithm 3.2 certifies and returns the current point x, not the trial xb.
 if native<=1e-6,status='native_pg_stop';break;end
 J=find(x~=0);
 omega=Inf;omegab=Inf;
 if ~isempty(J),omega=mu+min(-.25*P.alpha(J).*abs(x(J)).^(-1.5));end
 Jb=find(xb~=0);
 if ~isempty(Jb),omegab=mu+min(-.25*P.alpha(Jb).*abs(xb(Jb)).^(-1.5));end
 y=xb;branch='half_prox';shift=NaN;
 if isequal(sign(x),sign(xb))&&omega>=o.hpg_delta*omegab&&~isempty(J)
  gr=g(J)+.5*P.alpha(J).*sign(x(J))./sqrt(abs(x(J)));
  Hr=H(J,J)+diag(-.25*P.alpha(J).*abs(x(J)).^(-1.5));Hr=(Hr+Hr')/2;
  mineig=min(eig(Hr));shift=-o.hpg_b1*min(0,mineig)+o.hpg_b2*norm(gr)^o.hpg_sigma;
  q=zeros(P.d,1);q(J)=-(Hr+shift*eye(numel(J)))\gr;nlin=nlin+1;
  for j=0:o.external_backtracks-1
   t=.5^j;z=x+t*q;
   if np_difference(P,x,z)<=1e-4*t*(gr'*q(J)),y=z;branch='subspace_newton';break;end
  end
 end
 dn=norm(y-x);[~,Ry]=np_certificate(P,y,[],o);fy=np_eval(P,y);
 hist(end+1,:)=[k fy Ry dn NaN NaN NaN 0 numel(J) shift NaN NaN NaN 0];
 log{end+1}=struct('mu',mu,'native_pg_residual',native,'branch',branch,'shift',shift);
 x=y;X(:,end+1)=x;elapsed(end+1)=toc(timer);
 if dn==0,status='zero_step';break;end
end
out=np_external_finish(P,x,o,'HpgSRN',status,native,X,hist,elapsed,log,timer);
out.counts.linear_solves=nlin;out.variant='Algorithm_3_2_exact_reduced_solve';
out.parameters=struct('alpha',o.hpg_alpha,'delta',o.hpg_delta,'b1',o.hpg_b1, ...
 'b2',o.hpg_b2,'sigma',o.hpg_sigma,'armijo',1e-4,'beta',.5,'BB_bounds',[1e-20 1e20]);
end
