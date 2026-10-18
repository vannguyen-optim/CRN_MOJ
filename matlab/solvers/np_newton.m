function [q,ok,s]=np_newton(P,x,o,kind)
% Face linear direction first; support convex QP is the second option.
% Both use the true NEGATIVE fractional-power curvature on actual nonzeros.
if nargin<4,kind='support';end
assert(strcmp(o.curvature_modification,'absolute_eigenvalues'),'Unknown curvature modification');
[~,g,H,~,margin]=np_eval(P,x);g=g+np_power_grad(P,x);
q=zeros(P.d,1);s=struct('dimension',0,'shift',0,'mineig',NaN,'error',NaN, ...
 'qp',0,'linear',0,'kind',kind,'gradient_norm',NaN,'epsilon_target',NaN, ...
 'qp_reduced_gradient_norm',NaN,'reason','');
ok=false;
if strcmp(P.kind,'svm')&&any(margin==0)
 s.reason='unavailable_curvature_at_margin_kink';return;
end
Z=find(P.alpha>0&x==0);K=setdiff((1:P.d)',Z);
if strcmp(kind,'face'),[CC,~]=np_active_face(P,x);
elseif strcmp(kind,'support'),CC=P.C;
else,error('Unknown Newton direction kind: %s',kind);end
U0=np_null(CC(:,K),numel(K));U=zeros(P.d,size(U0,2));U(K,:)=U0;
s.dimension=size(U,2);ok=true;if isempty(U),s.reason='zero_dimensional';return;end
J=P.alpha>0&x~=0;D=zeros(P.d,1);
D(J)=P.alpha(J).*P.p(J).*(P.p(J)-1).*abs(x(J)).^(P.p(J)-2);
HR=U'*(H+diag(D))*U;HR=(HR+HR')/2;gr=U'*g;
if any(~isfinite(HR(:)))||any(~isfinite(gr))
 ok=false;s.reason='unavailable_numeric_curvature';return;
end
[VV,EE]=eig(HR);eigenvalues=diag(EE);s.mineig=min(eigenvalues);
modified=max(o.mN,abs(eigenvalues));
if all(eigenvalues>=o.mN),HH=HR;else,HH=VV*diag(modified)*VV';HH=(HH+HH')/2;end
s.shift=max(modified-eigenvalues);s.modification='absolute_eigenvalues';
s.true_eigenvalues=eigenvalues;s.modified_eigenvalues=modified;
s.gradient_norm=norm(gr);s.epsilon_target=min(.1,norm(gr));
if strcmp(kind,'face')||(isempty(P.A)&&all(isinf(P.lb))&&all(isinf(P.ub)))
 qr=-HH\gr;s.linear=1;
 % Finite-budget iterative refinement against the actual shifted system.
 for j=1:3
  if ~strcmp(kind,'face'),break;end
  err=HH*qr+gr;
  if norm(err)<=s.epsilon_target*norm(gr),break;end
  qr=qr-HH\err;s.linear=s.linear+1;
 end
 ok=all(isfinite(qr));
 if strcmp(kind,'face'),ok=ok&&norm(HH*qr+gr)<=s.epsilon_target*norm(gr);end
 if ~ok,s.reason='linear_system_accuracy_failure';end
else
 AR=[P.A*U;U(isfinite(P.ub),:);-U(isfinite(P.lb),:)];
 br=[P.b-P.A*x;P.ub(isfinite(P.ub))-x(isfinite(P.ub)); ...
     x(isfinite(P.lb))-P.lb(isfinite(P.lb))];
 [qr,ok,qmeta]=np_qp(HH,gr,AR,br,[],[],[],[],zeros(size(U,2),1),o);
 s.qp=double(qmeta.performed);
 if ~ok,s.reason=qmeta.message;end
end
if ok
 q=U*qr;
 if s.qp>0,s.qp_reduced_gradient_norm=norm(HH*qr+gr)/max(norm(gr),realmin);
 else,s.error=norm(HH*qr+gr)/max(norm(gr),realmin);end
end
end
