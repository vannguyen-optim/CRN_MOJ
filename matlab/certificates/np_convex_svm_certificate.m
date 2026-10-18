function [u,R,ok,work]=np_convex_svm_certificate(P,x,penalty,o)
% Minimum ORIGINAL CONVEX constrained subgradient norm, not a fractional certificate.
% l1 zero subgradients belong to [-1,1]; they are never free multipliers.
[~,q]=np_convex_svm_eval(P,x,penalty);x=x(:);work=0;
Z=[];
if strcmp(penalty,'l1')
 J=find(x(1:end-1)~=0);q(J)=q(J)+sign(x(J));Z=find(x(1:end-1)==0);
end
ok=np_feas(P,x)<=o.feas_tol;
if ~ok,u=nan(P.d,1);R=Inf;return;end
if isempty(P.C)&&isempty(P.A)
 u=q;u(Z)=u(Z)-min(1,max(-1,u(Z)));
 for j=1:P.d
  if P.lb(j)==P.ub(j),u(j)=0;
  elseif x(j)==P.lb(j),u(j)=min(0,u(j));
  elseif x(j)==P.ub(j),u(j)=max(0,u(j));end
 end
else
 U=np_null(P.C,P.d);qr=U'*q;E=eye(P.d);
 normals=P.A(P.A*x==P.b,:)';lo=find(x==P.lb&isfinite(P.lb));hi=find(x==P.ub&isfinite(P.ub));
 normals=[normals -E(:,lo) E(:,hi)];nn=size(normals,2);
 B=U'*[normals E(:,Z)];
 if isempty(B)||isempty(U),u=U*qr;
 else
  lb=[zeros(nn,1);-ones(numel(Z),1)];ub=[inf(nn,1);ones(numel(Z),1)];
  [v,solved,meta]=np_qp(B'*B,B'*qr,[],[],[],[],lb,ub,zeros(size(B,2),1),o);
  work=double(meta.performed);ok=ok&&solved;
  if isempty(v)||~all(isfinite(v))||~solved,u=nan(P.d,1);R=Inf;return;end
  v=min(ub,max(lb,v));u=U*(qr+B*v);
 end
end
R=norm(u);
end
