function [C,d]=np_qp_equalities(C,d)
if isempty(C),C=zeros(0,size(C,2));d=zeros(0,1);return;end
[U,S,V]=svd(full(C),'econ');ss=diag(S);r=sum(ss>max(size(C))*eps(max([ss;1])));
assert(norm(d-U(:,1:r)*(U(:,1:r)'*d),Inf)<=1e-9,'Inconsistent equalities');
d=(U(:,1:r)'*d)./ss(1:r);C=V(:,1:r)';
end
