function [x,ok,meta]=np_qp(H,f,A,b,C,d,lb,ub,x0,o)
% Native MathWorks QP. Solver flags do not certify original stationarity.
backend=np_native_optim();
meta=struct('flag',-99,'message','','identifier','','stack',[], ...
 'performed',false,'backend',backend.quadprog,'options_backend',backend.optimoptions);
x=x0;ok=false;
if np_work_budget('remaining')==0
 meta.message='Per-run QP solve budget exhausted';meta.flag=-98;return;
end
try
 opts=optimoptions('quadprog','Algorithm','active-set','Display','off', ...
  'MaxIterations',o.qp_max_iter,'ConstraintTolerance',o.qp_tol, ...
  'OptimalityTolerance',o.qp_tol,'StepTolerance',1e-14);
 HH=(H+H')/2;
 if ~np_work_budget('reserve')
  meta.message='Per-run QP solve budget exhausted';meta.flag=-98;return;
 end
 meta.performed=true;
 [x,~,flag,out]=quadprog(HH,f,A,b,C,d,lb,ub,x0,opts);
 meta.flag=flag;meta.output=out;
 if isfield(out,'message'),meta.message=out.message;end
 ok=flag>0&&~isempty(x)&&all(isfinite(x));
catch ME
 meta.identifier=ME.identifier;meta.message=ME.message;meta.stack=ME.stack;
end
end
