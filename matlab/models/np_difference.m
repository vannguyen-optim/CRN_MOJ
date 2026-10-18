function delta=np_difference(P,x,y)
% Cancellation-safe Phi(y)-Phi(x); evaluates the same original objective.
x=x(:);y=y(:);s=y-x;[~,g]=np_eval(P,x);
switch P.kind
 case 'lowdim'
  delta=g'*s+(P.lambda(1)+P.lambda(2))*sum(s.^2)+ ...
   P.lambda(3)*sum((6*x.^2-2).*s.^2+4*x.*s.^3+s.^4);
 case 'regression'
  delta=g'*s+.5*s'*P.B*s;
 case 'svm'
  margin=1-P.y.*(P.Xa*x);dm=-P.y.*(P.Xa*s);
  r=max(0,margin);rn=max(0,margin+dm);dr=rn-r;
  active=margin>0&margin+dm>0;dr(active)=dm(active);
  delta=sum(P.c.*(2*r.*dr+dr.^2));
 otherwise,error('Unknown problem kind');
end
J=find(P.alpha>0);ax=abs(x(J));ay=abs(y(J));p=P.p(J);
dp=ay.^p-ax.^p;ratio=zeros(size(ax));I=ax>0&ay>0;
ratio(I)=(ay(I)-ax(I))./ax(I);
% Large changes have no subtraction cancellation. In particular do not use
% log1p at a ratio rounded to -1: a tiny positive ay can have substantial ay^p.
local=I&abs(ratio)<.5;
dp(local)=ax(local).^p(local).*expm1(p(local).*log1p(ratio(local)));
delta=delta+sum(P.alpha(J).*dp);
end
