function out=np_solve(P,x,method,o)
% CRN, matched Newton ablation, and internal stabilized IRL1 diagnostic.
assert(any(strcmp(method,{'hybrid','no_newton','irl1'})), ...
 'Unknown method. External published solvers require a separate audited adapter.');
timer=tic;elapsed=0;np_work_budget('start',o.max_qp);x=x(:);hist=[];X=x;branches={};innerlog={};
triallog={};polishing={};face_history={};[~,face0]=np_active_face(P,x);face_history{1}=face0;
counts=struct('newton_attempts',0,'face_attempts',0,'support_attempts',0, ...
 'full',0,'damped',0,'face_full',0,'face_damped',0,'support_full',0,'support_damped',0, ...
 'direct_failed',0,'explore_accepted',0,'fallback_calls',0,'candidate_qp',0,'fallback_qp',0, ...
 'certificate_qp',0,'certificate_calls',1,'model_certificate_calls',0,'linear_solves',0, ...
 'model_steps',0,'rejected',0,'constructed_fallback',0,'refinement_failed',0);
status='outer_budget';[u,R,certok,cw]=np_certificate(P,x,[],o);counts.certificate_qp=cw;
initial_residual=R;f=np_eval(P,x);sigprev=[];newton=strcmp(method,'hybrid');initial_stationary=certok&&R<=o.tol;
if ~certok,status='initial_certificate_failure';end
for k=0:o.max_outer-1
 if ~certok,break;end
 % Explicit zero-start diagnostic: optional exploration before the first stop.
 diagnostic_first=k==0&&o.explore&&o.initial_explore_stationary&& ...
    any(strcmp(method,{'hybrid','no_newton'}));
 if R<=o.tol&&~diagnostic_first,status='success';break;end
 sig=sign(x(P.alpha>0));accepted=false;branch='';
 ns=struct('dimension',0,'shift',0,'mineig',NaN,'error',NaN);
 prox=NaN;inratio=NaN;y=x;uy=u;
 if newton&&k>=1&&isequal(sig,sigprev)
  counts.newton_attempts=counts.newton_attempts+1;
  kinds={'support'};if o.face_newton,kinds={'face','support'};end
  for ni=1:numel(kinds)
   kind=kinds{ni};counts.([kind '_attempts'])=counts.([kind '_attempts'])+1;
   [q,qok,ss]=np_newton(P,x,o,kind);
   counts.candidate_qp=counts.candidate_qp+ss.qp;counts.linear_solves=counts.linear_solves+ss.linear;
   if qok
    for j=0:o.newton_trials-1
     y=x+o.beta^j*q;raw_y=y;Z=find(P.alpha>0&x==0);y(Z)=0;
     % Restore literal active box bounds for the face direction only.
     if strcmp(kind,'face')
      lo=x==P.lb&isfinite(P.lb);hi=x==P.ub&isfinite(P.ub);y(lo)=P.lb(lo);y(hi)=P.ub(hi);
     end
     feasible=np_feas(P,y)<=o.feas_tol;J=P.alpha>0&x~=0;
     signs=all(sign(y(J))==sign(x(J)));yok=false;passed=false;
     if feasible&&signs
      [uy,~,yok,cw]=np_certificate(P,y,[],o);
      counts.certificate_qp=counts.certificate_qp+cw;counts.certificate_calls=counts.certificate_calls+1;
      passed=yok&&np_direct(P,x,y,uy,o);
     end
     triallog{end+1}=struct('iteration',k+1,'kind',kind,'trial',j, ...
      'step',o.beta^j,'feasible',feasible,'sign_preserving',signs,'certificate_ok',yok, ...
      'accepted',passed,'system',ss,'projection_norm',norm(y-raw_y));
     if any(y~=raw_y)
      triallog{end}.raw_output=raw_y;triallog{end}.polished_output=y;
     end
     if passed
      accepted=true;branch=[kind '_newton'];ns=ss;
      counts.full=counts.full+(j==0);counts.damped=counts.damped+(j>0);
      counts.([kind '_full'])=counts.([kind '_full'])+(j==0);
      counts.([kind '_damped'])=counts.([kind '_damped'])+(j>0);break;
     end
     counts.direct_failed=counts.direct_failed+1;counts.rejected=counts.rejected+1;
    end
   else
    triallog{end+1}=struct('iteration',k+1,'kind',kind,'trial',NaN, ...
     'step',NaN,'feasible',false,'sign_preserving',false,'certificate_ok',false, ...
     'accepted',false,'system',ss);
   end
   if accepted,break;end
  end
 end
 if any(strcmp(method,{'hybrid','no_newton'}))
  if ~accepted&&o.explore
   [y,uy,eok,ss]=np_explore(P,x,o);
   counts.candidate_qp=counts.candidate_qp+ss.qp;counts.certificate_qp=counts.certificate_qp+ss.certificate_qp;
   counts.certificate_calls=counts.certificate_calls+ss.certificate_calls;
   counts.model_certificate_calls=counts.model_certificate_calls+ss.model_certificate_calls;
   counts.model_steps=counts.model_steps+ss.steps;polishing=[polishing ss.polishing];
   if eok&&np_direct(P,x,y,uy,o)
    accepted=true;branch='exploration';counts.explore_accepted=counts.explore_accepted+1;
   else,counts.direct_failed=counts.direct_failed+1;counts.rejected=counts.rejected+1;end
  end
 else
  % A diagnostic IRL1 step, with original objective backtracking.
  T=P.L+o.eta+1;[~,g]=np_eval(P,x);
  w=P.alpha.*P.p.*max(abs(x),o.rw_eps).^(P.p-1);
  for tr=1:o.newton_trials
   [y,~,eok,mw,mcw,pp]=np_l1_model(P,x,g,T,w,[],o);
   counts.candidate_qp=counts.candidate_qp+mw;counts.certificate_qp=counts.certificate_qp+mcw;
   counts.model_certificate_calls=counts.model_certificate_calls+double(pp.model_certificate_called);
   counts.model_steps=counts.model_steps+1;if ~isempty(pp.raw),polishing{end+1}=pp;end
   if eok&&np_difference(P,x,y)<=0,accepted=true;branch='irl1';break;end
   T=2*T;counts.rejected=counts.rejected+1;
  end
  if accepted
   [uy,~,eok,cw]=np_certificate(P,y,[],o);
   counts.certificate_qp=counts.certificate_qp+cw;counts.certificate_calls=counts.certificate_calls+1;
   if ~eok,accepted=false;end
  end
 end
 if ~accepted&&strcmp(method,'irl1'),status='diagnostic_step_failure';break;end
 if ~accepted
  counts.fallback_calls=counts.fallback_calls+1;
  [y,uy,accepted,fs]=np_fallback(P,x,o);
  counts.fallback_qp=counts.fallback_qp+fs.model_qp;counts.model_steps=counts.model_steps+fs.inner;
  counts.model_certificate_calls=counts.model_certificate_calls+fs.model_certificate_calls;
  counts.certificate_calls=counts.certificate_calls+fs.certificate_calls;
  counts.certificate_qp=counts.certificate_qp+fs.certificate_qp;
  counts.constructed_fallback=counts.constructed_fallback+fs.constructed_used;
  counts.refinement_failed=counts.refinement_failed+fs.refinement_failed;
  polishing=[polishing fs.polishing];if o.inner_trace,innerlog{end+1}=fs;end
  if ~accepted,status=fs.reason;break;end
  branch='fallback';prox=fs.prox;inratio=fs.inner_ratio;
 end
 dn=norm(y-x);fy=np_eval(P,y);Ry=norm(uy);
 if dn>0,D=-np_difference(P,x,y)/dn^2;Q=Ry/dn;E=prox/dn;else,D=NaN;Q=NaN;E=NaN;end
 [~,facesig]=np_active_face(P,y);face_history{end+1}=facesig;
 hist=[hist;k+1 fy Ry dn D Q E 0 ns.dimension ns.shift ns.mineig ns.error inratio np_feas(P,y)];
 elapsed(end+1)=toc(timer);branches{end+1}=branch;sigprev=sig;x=y;u=uy;R=Ry;f=fy;
 if o.save_iterates,X(:,end+1)=x;end
 if R<=o.tol,status='success';break;end
 if dn==0,status='zero_step_nonstationary';break;end
end
if ~strcmp(status,'success')&&np_work_budget('remaining')==0,status='qp_budget';end
termination_reason=status;
% Independent COMMON minimum original residual at the returned point.
[uf,Rfinal,cok,cw]=np_certificate(P,x,[],o);
counts.certificate_qp=counts.certificate_qp+cw;counts.certificate_calls=counts.certificate_calls+1;
if ~cok
 if np_work_budget('remaining')==0,status='qp_budget';else,status='final_certificate_failure';end
elseif Rfinal<=o.tol,status='success';end
[objective,~,HH,F,margin]=np_eval(P,x);J=P.alpha>0&x~=0;dd=zeros(P.d,1);
dd(J)=P.alpha(J).*P.p(J).*(P.p(J)-1).*abs(x(J)).^(P.p(J)-2);
Z=find(P.alpha>0&x==0);K=setdiff((1:P.d)',Z);
UU=np_null(P.C(:,K),numel(K));[FF,~]=np_active_face(P,x);VV=np_null(FF(:,K),numel(K));
mineig=NaN;face_mineig=NaN;
curvature_available=all(isfinite(HH(:)))&&all(isfinite(dd))&&~(strcmp(P.kind,'svm')&&any(margin==0));
if curvature_available
 if ~isempty(UU)
  RR=UU'*(HH(K,K)+diag(dd(K)))*UU;if all(isfinite(RR(:))),mineig=min(eig(RR));end
 end
 if ~isempty(VV)
  RR=VV'*(HH(K,K)+diag(dd(K)))*VV;if all(isfinite(RR(:))),face_mineig=min(eig(RR));end
 end
end
changes=0;last=0;
if o.save_iterates
 S=sign(X(P.alpha>0,:));ch=find(any(diff(S,1,2)~=0,1));changes=numel(ch);if ~isempty(ch),last=ch(end);end
end
fc=[];for j=2:numel(face_history),if ~isequal(face_history{j},face_history{j-1}),fc(end+1)=j-1;end,end
lastface=0;if ~isempty(fc),lastface=fc(end);end
counts.total_qp=np_work_budget('used');
assert(counts.candidate_qp+counts.fallback_qp+counts.certificate_qp==counts.total_qp, ...
 'QP component accounting mismatch');
out=struct('method',method,'status',status,'success',strcmp(status,'success'), ...
 'termination_reason',termination_reason,'x',x,'u',uf,'residual',Rfinal,'feasibility',np_feas(P,x), ...
 'objective',objective,'F',F,'selected',nnz(P.alpha>0&x~=0), ...
 'initial_residual',initial_residual,'elapsed',elapsed,'iterations',size(hist,1),'history',hist,'branches',{branches}, ...
 'iterates',X,'inner',{innerlog},'counts',counts,'seconds',toc(timer), ...
 'newton_trials',{triallog},'raw_polishing',{polishing},'face_history',{face_history}, ...
 'initial_stationary',initial_stationary,'initial_stationary_exploration', ...
 initial_stationary&&o.initial_explore_stationary&&o.explore&&any(strcmp(method,{'hybrid','no_newton'})), ...
 'support_changes',changes,'last_support_change',last,'reduced_mineig',mineig, ...
 'face_mineig',face_mineig,'final_curvature_available',curvature_available, ...
 'face_dimension',size(VV,2),'face_changes',numel(fc),'last_face_change',lastface, ...
 'local_regular_condition','not_assessed','accurate_reference_limit',[], ...
 'min_margin',min([abs(margin);Inf]));
end
