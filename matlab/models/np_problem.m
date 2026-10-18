function P=np_problem(kind,data,lambda,constraint)
P.kind=kind;P.lambda=lambda(:);P.A=[];P.b=[];
switch kind
 case 'lowdim'
  P.d=2;P.p=.5*ones(2,1);P.alpha=(.10*lambda(1)+.15*lambda(2)+.05*lambda(3))*ones(2,1);
  P.L=2*(lambda(1)+lambda(2))+44*lambda(3);P.lb=-2*ones(2,1);P.ub=2*ones(2,1);
  P.C=zeros(0,2);P.rhs=zeros(0,1);P.A=zeros(0,2);P.b=zeros(0,1);return;
 case 'regression'
  P.d=size(data.A1,2);P.A1=data.A1;P.A2=data.A2;P.b1=data.b1(:);P.b2=data.b2(:);
  P.p=.5*ones(P.d,1);P.alpha=lambda(3)*ones(P.d,1);
  P.B=lambda(1)/size(P.A1,1)*(P.A1'*P.A1)+lambda(2)/size(P.A2,1)*(P.A2'*P.A2);
  P.h=lambda(1)/size(P.A1,1)*(P.A1'*P.b1)+lambda(2)/size(P.A2,1)*(P.A2'*P.b2);
  P.L=norm(P.B,2);
 case 'svm'
  P.d=size(data.X,2)+1;P.Xa=[data.X ones(size(data.X,1),1)];P.y=data.y(:);
  assert(any(P.y==1)&&any(P.y==-1),'Both classes required');
  P.p=.5*ones(P.d,1);P.alpha=[lambda(1)*ones(P.d-1,1);0];
  P.c=zeros(size(P.y));P.c(P.y==1)=lambda(2)/nnz(P.y==1);P.c(P.y==-1)=lambda(3)/nnz(P.y==-1);
  % Spectral bound sharper than the trace bound, still a valid global Lipschitz constant.
  P.L=2*norm(P.Xa'*bsxfun(@times,P.Xa,P.c),2);
 otherwise,error('Unknown problem');
end
P.lb=-inf(P.d,1);P.ub=inf(P.d,1);P.C=zeros(0,P.d);P.rhs=zeros(0,1);P.A=zeros(0,P.d);P.b=zeros(0,1);
if strcmp(constraint,'box')
 P.lb=-5*ones(P.d,1);P.ub=5*ones(P.d,1);
elseif strcmp(constraint,'affine')
 P.C=data.C;P.rhs=zeros(size(P.C,1),1);
elseif ~strcmp(constraint,'free'),error('Unknown constraint');end
end
