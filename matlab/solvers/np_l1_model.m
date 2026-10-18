function [y,e,ok,work,certwork,polish] = np_l1_model(P,v,grad,T,w,Z,o)
% min grad'*(y-v)+T/2*||y-v||^2+sum w_j |y_j| over Omega with Z locked.
work=0;certwork=0;Z=unique(Z(:));
polish=struct('raw',[],'polished',[],'proposed_zeros',[],'projection_norm',0,'accepted',false, ...
 'model_certificate_called',false,'model_backend','','model_solver',[]);
if isempty(P.C) && isempty(P.A)
 y=v;e=nan(P.d,1);ok=false;
 if ~np_work_budget('reserve'),return;end
 work=1;target=v-grad/T;
 polish.model_backend='soft_threshold';
 y=sign(target).*max(0,abs(target)-w/T);y=min(P.ub,max(P.lb,y));y(Z)=0;
 ok=np_feas(P,y)<=o.feas_tol;
else
 K=setdiff((1:P.d)',Z);J=find(w(K)>0);nk=numel(K);nt=numel(J);
 y=zeros(P.d,1);
 if nk==0,e=zeros(P.d,1);ok=np_feas(P,y)<=o.feas_tol;return;end
 % Remove redundant affine rows using the original SVD representation.
 CK=P.C(:,K);[UU,SS,~]=svd(full(CK),'econ');sing=diag(SS);
 rk=sum(sing>max(size(CK))*eps(max([sing;1])));
 if rk>0,CC=UU(:,1:rk)'*CK;dd=UU(:,1:rk)'*P.rhs;
 else,CC=zeros(0,nk);dd=zeros(0,1);end
 if norm(P.rhs-UU(:,1:rk)*(UU(:,1:rk)'*P.rhs),Inf)>o.feas_tol
  e=nan(P.d,1);ok=false;return;
 end
 fast=strcmp(P.kind,'regression')&&~isempty(P.C)&&isempty(P.A) ...
  &&all(P.lb==-Inf)&&all(P.ub==Inf)&&nt==nk&&all(isfinite(w(K)));
 if fast
  [sol,ok,qmeta]=np_affine_l1_fast(T,v(K),grad(K),w(K),CC,dd,CK,P.rhs,o);
  polish.model_backend=qmeta.strategy;polish.model_solver=qmeta;
 else
  I=eye(nk);E=I(J,:);
  H=blkdiag(T*eye(nk),zeros(nt)); f=[grad(K)-T*v(K);w(K(J))];
  A=[P.A(:,K) zeros(size(P.A,1),nt); E -eye(nt); -E -eye(nt)];
  b=[P.b;zeros(2*nt,1)];
  x0=[v(K);abs(v(K(J)))];
  [sol,ok,qmeta]=np_qp(H,f,A,b,[CC zeros(rk,nt)],dd, ...
   [P.lb(K);zeros(nt,1)],[P.ub(K);inf(nt,1)],x0,o);
  polish.model_backend='native_active_set';
 end
 work=double(qmeta.performed);
 if ~ok,e=nan(P.d,1);return;end
 y(K)=sol(1:nk);
 polish.raw=y;
 % A small coefficient only PROPOSES a new face. Feasibility is restored first.
 ZZ=union(Z,find(P.alpha>0 & abs(y)<=o.zero_snap));
 polish.proposed_zeros=ZZ;
 [yp,faceok]=np_face(P,y,ZZ,o);
 if faceok,y=yp;else,[y,ok]=np_face(P,y,Z,o);end
 polish.polished=y;polish.projection_norm=norm(y-polish.raw);
 polish.accepted=faceok;
end
% Compute an actual convex model subgradient, independently of solver flags.
polish.model_certificate_called=true;
[e,certok,cw]=np_model_certificate(P,y,grad+T*(y-v),w,Z,o);
certwork=cw;ok=ok&&certok;
end
