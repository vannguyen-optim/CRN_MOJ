function record=np_execute(job,cfg)
backend=np_native_optim();
load(fullfile(cfg.output,job.input),'D');o=cfg.options;timer=tic;cpu_timer=cputime;
try
 kind=job.experiment;if any(strcmp(kind,{'wdbc','svm'})),kind='svm';end
 P=np_problem(kind,D,job.lambda,job.geometry);x=D.x0;
 if strcmp(job.geometry,'affine'),[x,ok]=np_face(P,x,[],o);assert(ok,'Initial projection failed');end
 model_seconds=toc(timer);solver_timer=tic;cpu_timer=cputime;
 switch job.method
  case 'CRN',out=np_solve(P,x,'hybrid',o);
  case 'CRN0',out=np_solve(P,x,'no_newton',o);
  case 'IReNA',out=np_irena(P,x,o);
  case 'HpgSRN',out=np_hpgsrn(P,x,o);
  case 'SSQP',out=np_ssqp(P,x,o);
  case 'l1_SVM',out=np_l1_svm(P,x,o);
  case 'l2_SVM',out=np_l2_svm(P,x,o);
 end
 solver_seconds=toc(solver_timer);out.method=job.method;
 out.cpu_seconds=cputime-cpu_timer;out.seconds=model_seconds+solver_seconds;out.model_seconds=model_seconds;
 out.initial_point=x;out.success=out.success&&out.feasibility<=o.feas_tol;
 if strcmp(kind,'lowdim')
  ref=np_low_reference(job.lambda);out.reference=ref;out.reference_gap=out.objective-ref.value;
 elseif strcmp(kind,'regression')
  out.mse1=mean((D.At1*out.x-D.bt1).^2);out.mse2=mean((D.At2*out.x-D.bt2).^2);
 else
  % Predictive evaluation occurs outside recorded optimization time.
  out.validation=np_class_metrics(D.Xval,D.yval,out.x);
 end
 if cfg.local_diagnostics&&strcmp(job.method,'CRN')&&out.success
  out.local_diagnostics=np_local_diagnostics(P,out,o);
  out.local_regular_condition=out.local_diagnostics.normal_condition;
  if isfield(out.local_diagnostics,'reference_certified')&&out.local_diagnostics.reference_certified
   out.accurate_reference_limit=out.local_diagnostics.reference;
  end
 end
 record=struct('job',job,'out',out,'options',o,'completed_at',datestr(now,30),'optimization_backend',backend);
catch ME
 out=struct('method',job.method,'status','execution_error','termination_reason','execution_error', ...
  'success',false,'seconds',toc(timer),'message',ME.message,'identifier',ME.identifier, ...
  'residual',Inf,'feasibility',Inf,'objective',NaN,'selected',NaN,'iterations',NaN, ...
  'x',D.x0,'history',zeros(0,14),'counts',struct('total_qp',NaN,'candidate_qp',NaN, ...
  'fallback_qp',NaN,'certificate_qp',NaN,'linear_solves',NaN));
 record=struct('job',job,'out',out,'options',o,'completed_at',datestr(now,30),'optimization_backend',backend,'error_stack',ME.stack);
end
path=fullfile(cfg.output,'runs',[job.id '.mat']);temporary=[path '.partial'];
save(temporary,'record','-v7');movefile(temporary,path,'f');
end
