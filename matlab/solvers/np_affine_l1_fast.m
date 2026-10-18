function [sol,ok,meta]=np_affine_l1_fast(T,v,grad,w,C,d,originalC,originalRhs,o)
% One logical equality-constrained weighted-L1 model solve.
% C,d are the SAME independent-row equations used by np_l1_model;
% originalC contains only the unlocked columns of the original equations.
% Sparse IP supplies a dual seed. The exact soft-threshold dual refinement
% must certify the original equations and full epigraph KKT conditions.
% A failed fast attempt falls back to the original native active-set solver
% inside this same reserved logical model solve.
v=v(:);grad=grad(:);w=w(:);d=d(:);originalRhs=originalRhs(:);
n=numel(v);m=size(C,1);sol=[v;abs(v)];ok=false;
backend=np_native_optim();
meta=struct('performed',false,'strategy','not_performed','flag',-99, ...
 'backend',backend.quadprog,'options_backend',backend.optimoptions, ...
 'backend_qp_calls',0,'fallback_used',false,'ip_flag',NaN, ...
 'active_flag',NaN,'ip_output',[],'active_output',[],'dual',[], ...
 'message','','identifier','','stack',[], ...
 'fast_error_identifier','','fast_error_message','', ...
 'equality_multiplier_coordinates','independent_rows');
if np_work_budget('remaining')==0
 meta.flag=-98;meta.message='Per-run QP solve budget exhausted';return;
end
try
 activeOptions=optimoptions('quadprog','Algorithm','active-set','Display','off', ...
  'MaxIterations',o.qp_max_iter,'ConstraintTolerance',o.qp_tol, ...
  'OptimalityTolerance',o.qp_tol,'StepTolerance',1e-14);
catch ME
 meta.identifier=ME.identifier;meta.message=ME.message;meta.stack=ME.stack;return;
end
if ~np_work_budget('reserve')
 meta.flag=-98;meta.message='Per-run QP solve budget exhausted';return;
end
meta.performed=true;
try
 assert(isscalar(T)&&isfinite(T)&&T>0&&n>0&&numel(grad)==n ...
  &&numel(w)==n&&all(isfinite(w))&&all(w>0), ...
  'CRN:AffineModel','The fast model requires positive finite unlocked weights.');
 I=speye(n);Hs=blkdiag(T*I,sparse(n,n));As=[I -I;-I -I];
 bs=zeros(2*n,1);Cs=[sparse(C) sparse(m,n)];
 linear=[grad-T*v;w];lower=[-inf(n,1);zeros(n,1)];upper=inf(2*n,1);
 sparseOptions=optimoptions('quadprog','Algorithm','interior-point-convex', ...
  'Display','off','MaxIterations',o.qp_max_iter, ...
  'ConstraintTolerance',o.qp_tol,'OptimalityTolerance',o.qp_tol, ...
  'StepTolerance',1e-14,'LinearSolver','sparse');
 meta.backend_qp_calls=meta.backend_qp_calls+1;
 [seed,~,ipFlag,ipOutput,ipLambda]=quadprog(Hs,linear,As,bs,Cs,d, ...
  lower,upper,[],sparseOptions);
 meta.ip_flag=ipFlag;meta.ip_output=ipOutput;
 assert(~isempty(seed)&&all(isfinite(seed))&&isstruct(ipLambda) ...
  &&isfield(ipLambda,'eqlin')&&numel(ipLambda.eqlin)==m ...
  &&all(isfinite(ipLambda.eqlin)), ...
  'CRN:AffineSeed','The sparse solver did not return a finite equality-dual seed.');
 [candidate,~,dual]=dualRefineAffine(T,v,grad,w,C,d,originalC,originalRhs, ...
  ipLambda.eqlin,o.feas_tol,o.qp_tol);
 meta.dual=dual;
 if dual.converged&&numel(candidate)==2*n&&all(isfinite(candidate)) ...
   &&dual.reducedFeasibility<=dual.target ...
   &&dual.originalFeasibility<=dual.target&&dual.scaledKKT<=dual.target
  sol=candidate;ok=true;meta.strategy='sparse_ip_dual';meta.flag=1;
  meta.message='Original affine feasibility and epigraph KKT conditions certified';
  return;
 end
 meta.fast_error_identifier='CRN:AffineDualUncertified';
 meta.fast_error_message=dual.reason;
catch ME
 meta.fast_error_identifier=ME.identifier;meta.fast_error_message=ME.message;
end

% This direct call intentionally does not reserve a second logical QP.
% Its model, starting point, options and finite-output test match np_qp.
meta.fallback_used=true;meta.strategy='native_active_set_fallback';
try
 I=eye(n);H=blkdiag(T*I,zeros(n));A=[I -I;-I -I];
 linear=[grad-T*v;w];lower=[-inf(n,1);zeros(n,1)];upper=inf(2*n,1);
 meta.backend_qp_calls=meta.backend_qp_calls+1;
 [sol,~,flag,output]=quadprog((H+H')/2,linear,A,zeros(2*n,1), ...
  [C zeros(m,n)],d,lower,upper,[v;abs(v)],activeOptions);
 meta.flag=flag;meta.active_flag=flag;meta.active_output=output;
 if isfield(output,'message'),meta.message=output.message;end
 ok=flag>0&&~isempty(sol)&&all(isfinite(sol));
catch ME
 meta.identifier=ME.identifier;meta.message=ME.message;meta.stack=ME.stack;
end
end

function [z,lambda,meta] = dualRefineAffine(T,v,grad,w,C,d,originalC,originalRhs,nu0,feasTol,qpTol)
% Certified small-dual refinement used only by the affine model backend.
% Solve the SAME equality-constrained weighted-L1 proximal model by refining
% its small dual. No experiment source, work counter, tolerance, or file is
% changed here. All output certificates refer to the original epigraph QP.

v=v(:); grad=grad(:); w=w(:); d=d(:); originalRhs=originalRhs(:);
n=numel(v); m=size(C,1);
assert(isscalar(T)&&isfinite(T)&&T>0,'T must be positive and finite.');
assert(numel(grad)==n&&numel(w)==n&&size(C,2)==n&&numel(d)==m, ...
    'Inconsistent reduced-model dimensions.');
assert(size(originalC,2)==n&&numel(originalRhs)==size(originalC,1), ...
    'Inconsistent original affine dimensions.');
assert(all(isfinite(v))&&all(isfinite(grad))&&all(isfinite(w))&&all(w>=0) ...
    &&all(isfinite(C(:)))&&all(isfinite(d)) ...
    &&all(isfinite(originalC(:)))&&all(isfinite(originalRhs)), ...
    'The model data must be finite and the weights nonnegative.');
assert(isscalar(feasTol)&&isfinite(feasTol)&&feasTol>0 ...
    &&isscalar(qpTol)&&isfinite(qpTol)&&qpTol>0, ...
    'The supplied tolerances must be positive and finite.');
if isempty(nu0), nu=zeros(m,1); else, nu=nu0(:); end
assert(numel(nu)==m&&all(isfinite(nu)),'Invalid initial equality multiplier.');

a=v-grad/T;
threshold=w/T;
linearY=grad-T*v;
target=min([1e-12,0.1*qpTol,0.01*feasTol]);
maxIterations=200;
maxBacktracks=60;
armijo=1e-4;
meta=struct('converged',false,'iters',0,'iterations',0, ...
    'residual',Inf,'reducedFeasibility',Inf,'originalFeasibility',Inf, ...
    'scaledKKT',Inf,'target',target,'reason','iteration_limit', ...
    'newtonSteps',0,'gradientSteps',0,'lineSearchEvaluations',0, ...
    'lastRegularization',0,'dualValue',Inf,'nu',nu);

current=evaluate(nu);
best=current;
if acceptable(current)
    meta.converged=true; meta.reason='certified';
elseif m==0
    meta.reason='original_affine_equations_not_satisfied';
else
    Ldual=norm(full(C*C'),2)/T;
    if ~isfinite(Ldual)||Ldual<=0
        meta.reason='no_reduced_dual_direction';
    else
        for iteration=1:maxIterations
            % Generalized derivative: soft_0 is the identity, including at 0.
            active=(abs(current.shifted)>threshold)|(w==0);
            H=full(C(:,active)*C(:,active)')/T;
            H=(H+H')/2;
            [step,newtonOK,regularization]=newtonDirection(H,current.gradient);
            meta.lastRegularization=regularization;
            accepted=false;
            for directionKind=1:2
                if directionKind==1
                    if ~newtonOK, continue; end
                    direction=step;
                else
                    direction=-current.gradient/Ldual;
                end
                slope=current.gradient'*direction;
                if any(~isfinite(direction))||~isfinite(slope)||slope>=0
                    continue;
                end
                alpha=1;
                for backtrack=1:maxBacktracks
                    trial=evaluate(current.nu+alpha*direction);
                    meta.lineSearchEvaluations=meta.lineSearchEvaluations+1;
                    if trial.finite&&trial.residual<best.residual, best=trial; end
                    % A complete primal/KKT certificate takes precedence over
                    % a dual-merit comparison lost to roundoff near the root.
                    if acceptable(trial)
                        current=trial; best=trial; accepted=true;
                        meta.converged=true; meta.reason='certified';
                        break;
                    end
                    if trial.finite&&trial.value<=current.value+armijo*alpha*slope
                        if isequal(trial.nu,current.nu)
                            break;
                        end
                        current=trial; accepted=true;
                        break;
                    end
                    alpha=alpha/2;
                end
                if accepted
                    if directionKind==1, meta.newtonSteps=meta.newtonSteps+1;
                    else, meta.gradientSteps=meta.gradientSteps+1; end
                    break;
                end
            end
            if ~accepted
                meta.reason='dual_line_search_stalled';
                break;
            end
            meta.iters=iteration;
            meta.iterations=iteration;
            if current.residual<best.residual, best=current; end
            if meta.converged, break; end
        end
    end
end

if ~meta.converged, current=best; end
z=[current.y;abs(current.y)];
lambda=current.lambda;
meta.residual=current.residual;
meta.reducedFeasibility=current.reducedFeasibility;
meta.originalFeasibility=current.originalFeasibility;
meta.scaledKKT=current.scaledKKT;
meta.dualValue=current.value;
meta.nu=current.nu;

    function yes=acceptable(state)
        yes=state.finite&&state.reducedFeasibility<=target ...
            &&state.originalFeasibility<=target&&state.scaledKKT<=target;
    end

    function state=evaluate(multiplier)
        shifted=a-C'*multiplier/T;
        y=sign(shifted).*max(abs(shifted)-threshold,0);
        gradient=d-C*y;
        dualValue=T/2*(y'*y)+d'*multiplier;
        % Valid epigraph multipliers, with zero lower/upper-bound multipliers.
        % Tiny roundoff violations are measured below, never silently clipped.
        subgradient=T*(v-y)-grad-C'*multiplier;
        plus=(w+subgradient)/2;
        minus=(w-subgradient)/2;
        eqReduced=norm(C*y-d,Inf);
        eqOriginal=norm(originalC*y-originalRhs,Inf);
        epigraph=abs(y);
        stationY=T*y+linearY+plus-minus+C'*multiplier;
        stationT=w-plus-minus;
        stationarity=norm([stationY;stationT],Inf);
        dualViolation=max([0;-plus;-minus]);
        complementarity=max(abs([0;plus.*(y-epigraph);minus.*(-y-epigraph)]));
        scale=1+norm([linearY;w],Inf)+norm(T*y,Inf);
        scaledKKT=max([stationarity,dualViolation,complementarity])/scale;
        residual=max([eqReduced,eqOriginal,scaledKKT]);
        finite=all(isfinite(y))&&all(isfinite(gradient)) ...
            &&isfinite(dualValue)&&isfinite(residual);
        multipliers=struct('ineqlin',[plus;minus],'eqlin',multiplier, ...
            'lower',zeros(2*n,1),'upper',zeros(2*n,1));
        state=struct('nu',multiplier,'shifted',shifted,'y',y, ...
            'gradient',gradient,'value',dualValue,'lambda',multipliers, ...
            'reducedFeasibility',eqReduced,'originalFeasibility',eqOriginal, ...
            'scaledKKT',scaledKKT,'residual',residual,'finite',finite);
    end
end

function [direction,ok,regularization]=newtonDirection(H,gradient)
% The unregularized Newton system is used whenever its Cholesky is valid.
direction=zeros(size(gradient)); ok=false; regularization=0;
if isempty(H)||any(~isfinite(H(:)))||any(~isfinite(gradient)), return; end
I=eye(size(H,1));
matrixScale=max(1,norm(H,Inf));
[V,E]=eig(H); eigenvalues=diag(E);
rankTolerance=max(size(H))*eps(max([abs(eigenvalues);1]));
positive=eigenvalues>rankTolerance;
if sum(positive)<size(H,1)
    % Do not rely on chol alone to detect a numerically rank-deficient face.
    % Preserve its Newton components and regularize only the null component.
    components=V(:,positive)'*gradient;
    nullGradient=gradient-V(:,positive)*components;
    direction=-V(:,positive)*(components./eigenvalues(positive));
    if norm(nullGradient)>100*eps*max(1,norm(gradient))
        regularization=1e-8*matrixScale;
        direction=direction-nullGradient/regularization;
    end
    slope=gradient'*direction;
    if all(isfinite(direction))&&isfinite(slope)&&slope<0
        ok=true; return;
    end
end
for attempt=0:12
    [R,p]=chol(H+regularization*I);
    if p==0
        direction=-(R\(R'\gradient));
        slope=gradient'*direction;
        if all(isfinite(direction))&&isfinite(slope)&&slope<0
            ok=true; return;
        end
    end
    if regularization==0
        regularization=10*eps(matrixScale);
    else
        regularization=10*regularization;
    end
end
end
