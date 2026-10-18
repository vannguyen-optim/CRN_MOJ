function [e,ok,work]=np_model_certificate(P,y,g,w,Z,o)
% Actual convex-model subgradient; solver flags are never certificates.
work=0;ok=np_feas(P,y)<=o.feas_tol;e=g;
if ~ok,e(:)=NaN;return;end
K=setdiff((1:P.d)',Z);nz=K(y(K)~=0);zz=K(y(K)==0&w(K)>0);
e(nz)=e(nz)+w(nz).*sign(y(nz));
if isempty(P.C)&&isempty(P.A)
 e(Z)=0;
 for j=zz',e(j)=sign(e(j))*max(0,abs(e(j))-w(j));end
 for j=K'
  if P.lb(j)==P.ub(j),e(j)=0;
  elseif y(j)==P.lb(j),e(j)=min(0,e(j));
  elseif y(j)==P.ub(j),e(j)=max(0,e(j));end
 end
 return;
end
U0=np_null(P.C(:,K),numel(K));U=zeros(P.d,size(U0,2));U(K,:)=U0;
E=eye(P.d);N=P.A(P.A*y==P.b,:)';
lo=find(y==P.lb&isfinite(P.lb));hi=find(y==P.ub&isfinite(P.ub));
N=[N -E(:,lo) E(:,hi)];B=U'*[E(:,zz) N];er=U'*e;
if isempty(B)||isempty(U),e=U*er;return;end
nn=size(N,2);lb=[-w(zz);zeros(nn,1)];ub=[w(zz);inf(nn,1)];
[v,solved,qmeta]=np_qp(B'*B,B'*er,[],[],[],[],lb,ub,zeros(size(B,2),1),o);
work=double(qmeta.performed);ok=ok&&solved;
if solved,v=min(ub,max(lb,v));e=U*(er+B*v);else,e(:)=NaN;end
end
