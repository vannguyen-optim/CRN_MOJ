function out=np_convex_svm(P,x,penalty,o)
% Convex squared-hinge SVM, ORIGINAL convex certificate, predictive comparison only.
% The bias is unpenalized and C+ / C- match P.lambda exactly after division by lambda1.
assert(strcmp(P.kind,'svm'),'Convex SVM predictors apply only to SVM jobs');
assert(P.lambda(1)>0,'Positive regularization weight is required');
timer=tic;np_work_budget('start',o.max_qp);x=x(:);d=P.d;n=numel(P.y);dw=d-1;
nt=0;if strcmp(penalty,'l1'),nt=dw;elseif ~strcmp(penalty,'l2'),error('Unknown penalty');end
nv=d+n+nt;c=P.c/P.lambda(1);H=zeros(nv);f=zeros(nv,1);
H(d+(1:n),d+(1:n))=diag(2*c);
if nt,f(d+n+(1:nt))=1;else,H(1:dw,1:dw)=eye(dw);end
A=[-bsxfun(@times,P.Xa,P.y) -eye(n) zeros(n,nt)];b=-ones(n,1);
if nt
 W=[eye(dw) zeros(dw,1)];A=[A;W zeros(dw,n) -eye(dw);-W zeros(dw,n) -eye(dw)];
 b=[b;zeros(2*dw,1)];
end
A=[A;P.A zeros(size(P.A,1),n+nt)];b=[b;P.b];
C=[P.C zeros(size(P.C,1),n+nt)];rhs=P.rhs;
% Remove redundant equality rows while retaining the original affine set.
if ~isempty(C)
 [UC,SC,VC]=svd(full(C),'econ');ss=diag(SC);rr=sum(ss>max(size(C))*eps(max(ss)));
 defect=rhs-UC(:,1:rr)*(UC(:,1:rr)'*rhs);
 assert(norm(defect,Inf)<=o.feas_tol,'Inconsistent affine representation');
 C=VC(:,1:rr)';rhs=(UC(:,1:rr)'*rhs)./ss(1:rr);
end
lb=[P.lb;zeros(n+nt,1)];ub=[P.ub;inf(n+nt,1)];
warm=[x;max(0,1-P.y.*(P.Xa*x));abs(x(1:nt))];
counts=struct('newton_attempts',0,'full',0,'damped',0,'explore_accepted',0, ...
 'fallback_calls',0,'candidate_qp',0,'fallback_qp',0,'certificate_qp',0, ...
 'certificate_calls',0,'model_certificate_calls',0,'linear_solves',0,'model_steps',0, ...
 'rejected',0,'face_attempts',0,'support_attempts',0,'face_full',0,'face_damped',0, ...
 'support_full',0,'support_damped',0,'constructed_fallback',0,'refinement_failed',0);
qp_log={};polishing={};lockedZ=[];u=nan(d,1);R=Inf;certok=false;status='qp_failure';qp_iterations=0;
for attempt=1:3
 [v,qok,qmeta]=np_qp(H,f,A,b,C,rhs,lb,ub,warm,o);
 counts.candidate_qp=counts.candidate_qp+double(qmeta.performed);qp_log{end+1}=qmeta;
 if isfield(qmeta,'output')&&isfield(qmeta.output,'iterations'),qp_iterations=qp_iterations+qmeta.output.iterations;end
 if isempty(v)||~all(isfinite(v)),status='qp_failure';break;end
 y=v(1:d);raw=y;
 if strcmp(penalty,'l1')
  Z=find(abs(y(1:end-1))<=1e-8);[proposal,fok]=np_face(P,y,Z,o);
  if fok,y=proposal;end
 end
 [uy,Ry,cok,cw]=np_convex_svm_certificate(P,y,penalty,o);
 counts.certificate_calls=counts.certificate_calls+1;counts.certificate_qp=counts.certificate_qp+cw;
 polishing{end+1}=struct('attempt',attempt,'raw',raw,'proposed',y,'norm',norm(y-raw), ...
  'original_convex_residual',Ry,'certificate_ok',cok,'qp_flag',qmeta.flag);
 x=y;u=uy;R=Ry;certok=cok;
 if certok&&R<=o.tol,status='success';break;end
 if ~qok,status='qp_not_converged';else,status='original_convex_certificate_not_met';end
 % Finite optional polish for free l1 models: propose a strict-complementarity
 % zero face, resolve the ORIGINAL QP, and recheck the FULL convex certificate.
 free_geometry=isempty(P.C)&&isempty(P.A)&&all(~isfinite(P.lb))&&all(~isfinite(P.ub));
 if attempt>=3||~strcmp(penalty,'l1')||~free_geometry,break;end
 [~,gs]=np_convex_svm_eval(P,x,'l1');Z=find(abs(x(1:end-1))<=1e-3&abs(gs(1:end-1))<1-1e-4);
 Z=setdiff(Z,lockedZ);if isempty(Z),break;end
 lockedZ=union(lockedZ,Z);EZ=eye(d);C=[C;EZ(Z,:) zeros(numel(Z),n+nt)];rhs=[rhs;zeros(numel(Z),1)];
 xwarm=x;xwarm(lockedZ)=0;warm=[xwarm;max(0,1-P.y.*(P.Xa*xwarm));abs(xwarm(1:nt))];
end
termination_reason=status;
% A final independent original CONVEX certificate; np_certificate is deliberately unused.
[uf,Rfinal,cok,cw]=np_convex_svm_certificate(P,x,penalty,o);
counts.certificate_calls=counts.certificate_calls+1;counts.certificate_qp=counts.certificate_qp+cw;
if cok&&Rfinal<=o.tol,status='success';
elseif np_work_budget('remaining')==0,status='qp_budget';
elseif ~cok,status='final_convex_certificate_failure';end
[value,~,regularizer]=np_convex_svm_eval(P,x,penalty);fractional=np_eval(P,x);
counts.total_qp=np_work_budget('used');
assert(counts.candidate_qp+counts.certificate_qp==counts.total_qp,'Convex QP accounting mismatch');
out=struct('method',[penalty '_svm'],'status',status,'common_status',status,'success',strcmp(status,'success'), ...
 'termination_reason',termination_reason,'x',x,'u',uf,'residual',Rfinal,'feasibility',np_feas(P,x), ...
 'objective',value,'F',[],'selected',nnz(x(1:end-1)~=0),'regularizer',regularizer, ...
 'iterations',NaN,'qp_iterations',qp_iterations,'history',zeros(0,14),'branches',{{}}, ...
 'iterates',x,'inner',{{}},'counts',counts,'seconds',toc(timer),'newton_trials',{{}}, ...
 'raw_polishing',{polishing},'qp_log',{qp_log},'support_changes',NaN,'last_support_change',NaN, ...
 'reduced_mineig',NaN,'min_margin',min(abs(1-P.y.*(P.Xa*x))), ...
 'model_kind',['convex_' penalty '_predictive'],'certificate_kind','original_convex_constrained', ...
 'scalarized_fractional_objective',fractional,'class_coefficients', ...
 [P.lambda(2)/(P.lambda(1)*nnz(P.y==1)),P.lambda(3)/(P.lambda(1)*nnz(P.y==-1))]);
end
