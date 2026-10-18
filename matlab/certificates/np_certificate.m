function [u,R,ok,work]=np_certificate(P,x,shift,o)
% Minimum-norm structured original certificate, in floating point.
% Only literal original active inequalities and penalized zeros are used.
% Eliminate free equality/zero multipliers before the nonnegative QP.
if nargin<3||isempty(shift),shift=zeros(P.d,1);end
[~,g]=np_eval(P,x);q=g+np_power_grad(P,x)+shift;
Z=find(P.alpha>0&x==0);K=setdiff((1:P.d)',Z);
work=0;ok=np_feas(P,x)<=o.feas_tol;
if ~ok,u=nan(P.d,1);R=Inf;return;end
if isempty(P.C)&&isempty(P.A)
 r=q;r(Z)=0;
 for j=K'
  if P.lb(j)==P.ub(j),r(j)=0;
  elseif x(j)==P.lb(j),r(j)=min(0,q(j));
  elseif x(j)==P.ub(j),r(j)=max(0,q(j));end
 end
else
 U0=np_null(P.C(:,K),numel(K));U=zeros(P.d,size(U0,2));U(K,:)=U0;
 qr=U'*q;E=eye(P.d);normals=P.A(P.A*x==P.b,:)';
 lo=find(x==P.lb&isfinite(P.lb));hi=find(x==P.ub&isfinite(P.ub));
 normals=[normals -E(:,lo) E(:,hi)];B=U'*normals;
 if isempty(B)||isempty(U),r=U*qr;
 else
  [v,solved,qmeta]=np_qp(B'*B,B'*qr,[],[],[],[],zeros(size(B,2),1), ...
      inf(size(B,2),1),zeros(size(B,2),1),o);
  work=double(qmeta.performed);ok=ok&&solved;
  if ~solved,u=nan(P.d,1);R=Inf;return;end
  v=max(0,v);r=U*(qr+B*v);
 end
end
u=r-shift;R=norm(r);
end
