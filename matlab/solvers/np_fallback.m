function [y,u,ok,s]=np_fallback(P,x,o)
% Tangent weights; locked zeros; strict relative model-error check.
% Manuscript (van_115) always constructs the proximal certificate.
% Proximal minimum-norm refinement is optional and budgeted.
[fx,~]=np_eval(P,x);v=x;T=P.L+o.eta+1;c=T-(P.L+o.eta)/2-o.chi;
assert(T>(P.L+o.eta)/2&&o.chi>=0&&c>0,'Invalid fallback parameters');
s=struct('model_qp',0,'certificate_qp',0,'inner',0,'zero_monotone',true, ...
 'reason','inner_budget','trace',[],'prox',Inf,'inner_ratio',Inf, ...
 'model_certificate_calls',0,'certificate_calls',0,'constructed_used',0, ...
 'refinement_failed',0,'polishing',{{}});
y=x;u=nan(P.d,1);ok=false;
for r=1:o.max_inner
 Z=find(P.alpha>0&v==0);J=P.alpha>0&v~=0;w=zeros(P.d,1);
 w(J)=P.alpha(J).*P.p(J).*abs(v(J)).^(P.p(J)-1);
 [fv,g0]=np_eval(P,v);g=g0+o.eta*(v-x);
 [yn,e,modelok,mw,mcw,polish]=np_l1_model(P,v,g,T,w,Z,o);
 s.model_qp=s.model_qp+mw;s.certificate_qp=s.certificate_qp+mcw;s.inner=r;
 s.model_certificate_calls=s.model_certificate_calls+double(polish.model_certificate_called);
 if ~isempty(polish.raw),s.polishing{end+1}=polish;end
 if ~modelok
  if np_work_budget('remaining')==0,s.reason='qp_budget';
  else,s.reason='inner_QP_or_certificate_failure';end
  return;
 end
 ds=norm(yn-v);en=norm(e);[fy,gy]=np_eval(P,yn);
 relative_ok=en<=o.chi*ds;
 step=yn-v;
 descent_ok=np_difference(P,v,yn)+o.eta*(v-x)'*step+(o.eta/2+c)*ds^2<=0;
 % Do not relax the relative test or objective test near numerical cancellation.
 if ~relative_ok||~descent_ok,s.reason='inner_relative_or_descent_failure';return;end
 s.zero_monotone=s.zero_monotone&&all(yn(Z)==0);
 % Add the exact change of power derivatives at surviving nonzero indices.
 delta=zeros(P.d,1);Jy=P.alpha>0&yn~=0;
 delta(Jy)=P.alpha(Jy).*P.p(Jy).*sign(yn(Jy)).* ...
     (abs(yn(Jy)).^(P.p(Jy)-1)-abs(v(Jy)).^(P.p(Jy)-1));
 constructed=e+(gy-g0)+o.eta*(yn-v)-T*(yn-v)+delta;
 shift=o.eta*(yn-x);un=constructed-shift;prox=norm(constructed);
 chosen='constructed';s.constructed_used=s.constructed_used+1;
 if o.proximal_refine&&np_work_budget('remaining')>0
  [ur,pr,certok,cw]=np_certificate(P,yn,shift,o);
  s.certificate_qp=s.certificate_qp+cw;s.certificate_calls=s.certificate_calls+1;
  if certok&&pr<=prox,un=ur;prox=pr;chosen='minimum_norm';s.constructed_used=s.constructed_used-1;
  elseif ~certok,s.refinement_failed=s.refinement_failed+1;end
 end
 s.prox=prox;s.inner_ratio=en/max(ds,realmin);
 if o.inner_trace
  s.trace=[s.trace;r ds en prox numel(Z) nnz(P.alpha>0&yn==0) strcmp(chosen,'constructed')];
 end
 if all(isfinite(un))&&prox<=o.rho*norm(yn-x)&&np_difference(P,x,yn)+o.eta/2*norm(yn-x)^2<=0
  y=yn;u=un;ok=true;s.reason='accepted';return;
 end
 v=yn;
end
end
