function out=np_irena(P,x,o)
% IReNA Algorithm 3, full residual subspaces, exact reduced linear solve.
% https://arxiv.org/abs/2407.17216v4
assert(strcmp(P.kind,'regression')&&isempty(P.C)&&all(isinf(P.lb)));
timer=tic;np_work_budget('start',o.max_qp);epsv=ones(P.d,1);X=x;
hist=[];elapsed=0;log={};xp=[];gp=[];status='iteration_budget';native=Inf;nlin=0;
for k=1:o.max_outer
 [~,R,cok]=np_certificate(P,x,[],o);if cok&&R<=o.tol,status='common_stop';break;end
 [~,g,H]=np_eval(P,x);w=.5*P.alpha./sqrt(abs(x)+epsv);
 psi=zeros(P.d,1);Z=x==0;psi(Z)=sign(g(Z)).*max(0,abs(g(Z))-w(Z));
 phi=g+w.*sign(x);phi(Z)=0;
 I=x>0&g+w>0;phi(I)=min(g(I)+w(I),max(x(I),g(I)-w(I)));
 I=x<0&g-w<0;phi(I)=max(g(I)-w(I),min(x(I),g(I)+w(I)));
 native=max(norm(phi),norm(psi));
 if native<=1e-6&&all(epsv(x~=0)<=1e-6),status='native_residual_stop';break;end
 choosePsi=norm(psi)>=o.irena_eta*norm(phi);
 if choosePsi,W=find(psi~=0);else,W=find(phi~=0);end
 if isempty(W)
  W=find(x~=0);choosePsi=false;
  if isempty(W),status='native_stationary';break;end
 end
 beta=1/max(P.L,1e-12);
 if ~isempty(xp)
  s=x-xp;r=g-gp;if s'*r>0,beta=min(1e12,max(1e-12,(s'*s)/(s'*r)));end
 end
 xp=x;gp=g;accepted=false;G0=smooth(P,x)+sum(w.*abs(x));
 for j=1:o.external_backtracks
  ys=x;v=x(W)-beta*g(W);ys(W)=sign(v).*max(0,abs(v)-beta*w(W));
  if smooth(P,ys)+sum(w.*abs(ys))<=G0-o.irena_xi/2*norm(ys-x)^2,accepted=true;break;end
  beta=.5*beta;
 end
 if ~accepted,status='ist_backtracking_limit';break;end
 y=ys;branch='IST';shift=NaN;acceptedQP=false;
 if isequal(sign(ys),sign(x))&&any(x(W)~=0)
  gr=g(W)+w(W).*sign(x(W));D=-.25*P.alpha(W).*(abs(x(W))+epsv(W)).^(-1.5);
  Hr=H(W,W)+diag(D);Hr=(Hr+Hr')/2;
  % The nominal shift uses the entire current support, as in the paper.
  nominal=1e-8+1e-4*norm(gr)^.5+min(.25*P.alpha(x~=0).*(abs(x(x~=0))+epsv(x~=0)).^(-1.5));
  shift=max(nominal,1e-8-min(eig(Hr)));q=-(Hr+shift*eye(numel(W)))\gr;nlin=nlin+1;
  F0=smooth(P,x)+sum(P.alpha.*sqrt(abs(x)+epsv));
  t=1;projected_once=false;
  for j=1:o.external_backtracks
   yq=x;yq(W)=sign(x(W)).*max(0,sign(x(W)).*(x(W)+t*q));
   if isequal(sign(yq),sign(x))&&~projected_once
    B=x(W).*q<0;if any(B),t=min(1,min(-x(W(B))./q(B)));end
    yq(W)=sign(x(W)).*max(0,sign(x(W)).*(x(W)+t*q));projected_once=true;
   end
   Fq=smooth(P,yq)+sum(P.alpha.*sqrt(abs(yq)+epsv));
   if Fq<=F0-o.irena_xi/2*norm(yq-x)^2
    deltaQP=F0-Fq;
    deltaIST=G0-(smooth(P,ys)+sum(w.*abs(ys)));
    if isequal(sign(yq),sign(x))||deltaQP>=o.irena_nu*deltaIST
     y=yq;branch='QP';acceptedQP=true;
    end
    break;
   end
   t=.5*t;
  end
 end
 I=y~=0;
 if acceptedQP,epsv(I)=min(.9*epsv(I),epsv(I).^2);
 else
  epsv=max(epsv,1e-8);
  if choosePsi,epsv(I)=.9*epsv(I);else,epsv(I)=.9*epsv(I).^1.1;end
 end
 epsv=max(realmin,epsv);
 dn=norm(y-x);[~,Ry]=np_certificate(P,y,[],o);fy=np_eval(P,y);
 hist(end+1,:)=[k fy Ry dn NaN NaN NaN 0 numel(W) shift NaN NaN NaN 0];
 log{end+1}=struct('branch',branch,'epsilon',epsv,'native_residual',native, ...
 'beta',beta,'shift',shift,'subspace',W,'chosen_psi',choosePsi);
 x=y;X(:,end+1)=x;elapsed(end+1)=toc(timer);
 if dn==0&&all(epsv(x~=0)<=1e-6),status='zero_step';break;end
end
out=np_external_finish(P,x,o,'IReNA',status,native,X,hist,elapsed,log,timer);
out.counts.linear_solves=nlin;out.variant='Algorithm_3_full_subspaces_exact_reduced_solve';
out.parameters=struct('eta',o.irena_eta,'eta_phi',1,'eta_psi',1,'gamma',.5, ...
 'xi',o.irena_xi,'nu',o.irena_nu,'epsilon_initial',1,'epsilon_floor_IST',1e-8);
end
function f=smooth(P,x)
% Evaluate the smooth loss directly rather than subtracting the penalty.
f=P.lambda(1)*norm(P.A1*x-P.b1)^2/(2*size(P.A1,1)) ...
 +P.lambda(2)*norm(P.A2*x-P.b2)^2/(2*size(P.A2,1));
end
