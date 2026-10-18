function [y,ok]=np_face(P,x,Z,o)
% Set proposed zero coordinates EXACTLY to zero, then solve affine equations
% using only the remaining columns. No free-zero certificate at small nonzeros.
y=x(:); Z=unique(Z(:));y(Z)=0; K=setdiff((1:P.d)',Z);
if ~isempty(P.C)
 if isempty(K)
  ok=norm(P.rhs,Inf)<=o.feas_tol;return;
 end
 CK=P.C(:,K); defect=P.rhs-CK*y(K);
 y(K)=y(K)+pinv(full(CK))*defect;
end
if isempty(P.C) && isempty(P.A)
 y=min(P.ub,max(P.lb,y)); y(Z)=0;
end
ok=np_feas(P,y)<=o.feas_tol;
end
