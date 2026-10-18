function d=np_local_diagnostics(P,out,o)
np_native_optim();
% Post-run analysis; separate from optimization timing/work counters.
x=out.x;Z=find(P.alpha>0&x==0);K=setdiff((1:P.d)',Z);
[A,face]=np_active_face(P,x);U0=np_null(A(:,K),numel(K));U=zeros(P.d,size(U0,2));U(K,:)=U0;
[~,g,H,~,margin]=np_eval(P,x);b=g+np_power_grad(P,x);
D=zeros(P.d,1);J=P.alpha>0&x~=0;D(J)=-.25*P.alpha(J).*abs(x(J)).^(-1.5);
HF=U'*(H+diag(D))*U;HF=(HF+HF')/2;mu=NaN;
if ~isempty(HF),mu=min(eig(HF));end
d=struct('face_dimension',size(U,2),'face_mineig',mu,'curvature_threshold',max(o.mN,2*o.aN), ...
 'curvature_condition',~isempty(HF)&&mu>max(o.mN,2*o.aN), ...
 'coefficient_tolerance',1e-8,'representation_tolerance',o.tol,'face',face, ...
 'normal_lower_bound',NaN,'representation_residual',NaN,'normal_condition',false, ...
 'reference',[],'reference_residual',Inf,'iterate_errors',[],'error_ratios',[], ...
 'quadratic_error_ratios',[],'seconds',NaN);
timer=tic;lo=x==P.lb&isfinite(P.lb);hi=x==P.ub&isfinite(P.ub);
N=[P.A(P.A*x==P.b,:);-eye_rows(P.d,lo);eye_rows(P.d,hi)];
E=[P.C;eye_rows(P.d,ismember((1:P.d)',Z))];V=np_null(E,P.d);
if isempty(N)
 d.representation_residual=norm(V'*b);d.normal_condition=d.representation_residual<=o.tol;
 d.normal_lower_bound=Inf;
else
 n=size(N,1);AA=V'*N';bb=-V'*b;
 opts=optimoptions('linprog','Display','off','ConstraintTolerance',1e-11,'OptimalityTolerance',1e-11);
 % Upper bound one handles redundant/unbounded positive representations.
 [s,~,flag]=linprog([zeros(n,1);-1],[-eye(n) ones(n,1)],zeros(n,1), ...
  [AA zeros(size(AA,1),1)],bb,[zeros(n,1);0],[inf(n,1);1],opts);
 if flag>0
  d.normal_lower_bound=s(end);d.representation_residual=norm(AA*s(1:n)-bb);
  d.normal_condition=s(end)>d.coefficient_tolerance&&d.representation_residual<=o.tol;
 end
end
if ~isempty(U)&&d.curvature_condition&&all(abs(margin)>0)
 % Refine on the final fixed support/face without altering the returned solver output.
 y=x;
 for k=1:50
  [~,gy,Hy]=np_eval(P,y);gy=gy+np_power_grad(P,y);
  dy=zeros(P.d,1);I=P.alpha>0&y~=0;dy(I)=-.25*P.alpha(I).*abs(y(I)).^(-1.5);
  gr=U'*gy;HH=U'*(Hy+diag(dy))*U;
  if norm(gr)<=1e-12,break;end
  if min(eig((HH+HH')/2))<=0,break;end
  q=-U*(HH\gr);accepted=false;
  for j=0:30
   z=y+.5^j*q;z(Z)=0;
   if all(sign(z(J))==sign(x(J)))&&np_feas(P,z)<=o.feas_tol&&np_difference(P,y,z)<=0
    accepted=true;y=z;break;
   end
  end
  if ~accepted,break;end
 end
 np_work_budget('start',Inf);[~,rr,ok]=np_certificate(P,y,[],o);
 d.reference=y;d.reference_residual=rr;d.reference_certified=ok&&rr<=1e-10;
 if d.reference_certified&&isfield(out,'iterates')
  err=sqrt(sum(bsxfun(@minus,out.iterates,y).^2,1));d.iterate_errors=err;
  d.error_ratios=err(2:end)./err(1:end-1);d.quadratic_error_ratios=err(2:end)./err(1:end-1).^2;
 end
elseif isempty(U)
 d.reference=x;d.reference_residual=out.residual;d.reference_certified=out.success;
end
d.seconds=toc(timer);
end
function E=eye_rows(n,I)
E=eye(n);E=E(I,:);
end
