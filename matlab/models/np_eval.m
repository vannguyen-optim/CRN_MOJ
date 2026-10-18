function [val,g,H,F,margin] = np_eval(P,x)
% Original objective; no epsilon smoothing is used here.
x=x(:); margin=[];
switch P.kind
 case 'lowdim'
  l=P.lambda; t=sum(abs(x).^P.p);
  F=[(x(1)-1)^2+(x(2)+1)^2+.10*t; ...
     (x(1)+1)^2+(x(2)-1)^2+.15*t; sum((x.^2-1).^2)+.05*t];
  g=2*l(1)*(x-[1;-1])+2*l(2)*(x-[-1;1])+4*l(3)*x.*(x.^2-1);
  H=diag(2*(l(1)+l(2))+l(3)*(12*x.^2-4));
 case 'regression'
  F=[norm(P.A1*x-P.b1)^2/(2*size(P.A1,1)); ...
     norm(P.A2*x-P.b2)^2/(2*size(P.A2,1));sum(abs(x).^P.p)];
  g=P.B*x-P.h; H=P.B;
 case 'svm'
  margin=1-P.y.*(P.Xa*x);r=max(0,margin);
  F=[sum(abs(x(1:end-1)).^P.p(1:end-1)); ...
     mean(r(P.y==1).^2);mean(r(P.y==-1).^2)];
  g=-2*P.Xa'*(P.c.*P.y.*r);
  if nargout>=3
   I=margin>0; H=2*P.Xa(I,:)'*bsxfun(@times,P.Xa(I,:),P.c(I));
  else,H=[];end
 otherwise, error('Unknown problem kind');
end
val=P.lambda(:)'*F;
end
