function out=np_ssqp(P,x,o)
% SSQP full constrained convex-QP variant, equation (5.11).
% https://arxiv.org/abs/1407.7205
% |x_j|^.5 is [alpha_j^2*x_j]_+^.5+[-alpha_j^2*x_j]_+^.5.
timer=tic;np_work_budget('start',o.max_qp);X=x;hist=[];elapsed=0;native_log={};
stages=floor(log(o.ssqp_epsilon)/log(o.ssqp_sigma));
mu=o.ssqp_epsilon/o.ssqp_sigma^stages;status='iteration_budget';native=Inf;candidate=0;certwork=0;
J=find(P.alpha>0);E=eye(P.d);As=[diag(P.alpha(J).^2)*E(J,:);-diag(P.alpha(J).^2)*E(J,:)];
[C,d]=np_qp_equalities(P.C,P.rhs);projection=pinv(C);
for k=1:o.max_outer
 [~,R,cok,cw]=np_certificate(P,x,[],o);certwork=certwork+cw;
 if cok&&R<=o.tol,status='common_stop';break;end
 [fx,g]=np_eval(P,x);loss=fx-sum(P.alpha.*sqrt(abs(x)));
 [reg,rg,B]=np_ssqp_smoothing(P,x,mu);gg=g+rg;
 z=x-gg;if isempty(C),z=min(P.ub,max(P.lb,z));else,z=z+projection*(d-C*z);end
 native=norm(z-x);
 while native<=mu&&mu>o.ssqp_epsilon
  mu=max(o.ssqp_epsilon,o.ssqp_sigma*mu);[reg,rg,B]=np_ssqp_smoothing(P,x,mu);gg=g+rg;
  z=x-gg;if isempty(C),z=min(P.ub,max(P.lb,z));else,z=z+projection*(d-C*z);end
  native=norm(z-x);
 end
 if native<=mu&&mu<=o.ssqp_epsilon,status='native_mu_stop';break;end
 if np_work_budget('remaining')<2,status='qp_budget';break;end
 t=As*x;neg=t < -mu;pos=t > 2*mu;
 A=[P.A;As(neg,:);-As(pos,:)];b=[P.b-P.A*x;mu*ones(nnz(neg),1);t(pos)/2];
 [q,qok,qmeta]=np_qp(diag(B+P.L),gg,A,b,C,zeros(size(d)), ...
  P.lb-x,P.ub-x,zeros(P.d,1),o);candidate=candidate+double(qmeta.performed);
 if ~qok,status='qp_failure';break;end
 y=x+q;[fy,gy]=np_eval(P,y);ly=fy-sum(P.alpha.*sqrt(abs(y)));
 regy=np_ssqp_smoothing(P,y,mu);modelgap=loss+g'*q+.5*P.L*(q'*q)-ly;
 decrease=loss+reg-ly-regy;
 maxrow=max([P.alpha.^2;0]);J0=max(8*sum(P.alpha.^4)+4*P.L,2*maxrow+2);
 guarantee=mu^3.5/J0;
 if modelgap < -2e-12*max([1 abs(loss) abs(ly)])||decrease-guarantee < -2e-12*max([1 abs(loss+reg) abs(ly+regy)])
  status='majorization_failure';break;
 end
 dn=norm(q);[~,Ry,cok,cw]=np_certificate(P,y,[],o);certwork=certwork+cw;
 hist(end+1,:)=[k fy Ry dn NaN NaN NaN 0 0 0 NaN NaN NaN np_feas(P,y)];
 native_log{end+1}=struct('mu',mu,'native_residual',native,'model_gap',modelgap, ...
 'smoothed_decrease',decrease,'guaranteed_decrease',guarantee,'qp',qmeta);x=y;X(:,end+1)=x;elapsed(end+1)=toc(timer);
 if dn==0,status='zero_step';break;end
end
out=np_external_finish(P,x,o,'SSQP',status,native,X,hist,elapsed,native_log,timer);
out.counts.candidate_qp=candidate;out.counts.certificate_qp=certwork+out.counts.certificate_qp;
assert(out.counts.total_qp==candidate+out.counts.certificate_qp,'SSQP QP accounting');
out.variant='SSQP_full_QP_equation_5_11';out.mu=mu;
out.parameters=struct('sigma',o.ssqp_sigma,'epsilon',o.ssqp_epsilon,'known_L',P.L);
end
