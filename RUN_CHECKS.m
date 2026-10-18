function checks=RUN_CHECKS()
% Deterministic mathematical, data, solver, and acceptance checks.
root=CRN_SETUP();
assert(exist('quadprog','file')==2,'Optimization Toolbox is required.');
cfg=CRN_CONFIG_PUBLICATION();o=cfg.options;o.max_outer=100;o.max_inner=300;o.max_qp=500;
checks={};np_work_budget('start',Inf);
P=np_problem('lowdim',struct(),[.3 .3 .4],'box');
[~,R,ok]=np_certificate(P,zeros(2,1),[],o);assert(ok&&R==0);pass('Original certificate at literal zero');
x=[.7;-1.2];h=1e-6;[~,g,H]=np_eval(P,x);g=g+np_power_grad(P,x);
fd=zeros(2,1);for j=1:2,e=zeros(2,1);e(j)=h;fd(j)=(np_eval(P,x+e)-np_eval(P,x-e))/(2*h);end
assert(norm(fd-g)<1e-7);pass('Original gradient');
HH=H+diag(-.25*P.alpha.*abs(x).^(-1.5));fdH=zeros(2);
for j=1:2,e=zeros(2,1);e(j)=h;[~,gp]=np_eval(P,x+e);[~,gm]=np_eval(P,x-e);fdH(:,j)=(gp+np_power_grad(P,x+e)-gm-np_power_grad(P,x-e))/(2*h);end
assert(norm(fdH-HH)<1e-6);pass('True negative penalty curvature');
[q,qok,qs]=np_newton(P,[.2;.2],o,'face');
assert(qok&&all(qs.modified_eigenvalues>=o.mN)&&qs.mineig<0&&all(isfinite(q)));
pass('Spectral Newton modification and forcing accuracy');
y=x+[1e-5;-2e-5];assert(abs(np_difference(P,x,y)-(np_eval(P,y)-np_eval(P,x)))<1e-12);pass('Original objective difference');
ref=np_low_reference([.05 .90 .05]);Pr=np_problem('lowdim',struct(),[.05 .90 .05],'box');
assert(abs(ref.value-np_eval(Pr,ref.point))<1e-12&&ref.value<=np_eval(Pr,[-1;1]));pass('Polynomial reference minima');
v=[-3;-.3;0;.4;3];a=.2;xprox=np_half_prox(v,a);
for j=1:numel(v)
 F=@(t).5*(t-v(j)).^2+a*sqrt(abs(t));grid=linspace(-4,4,20001);
 assert(F(xprox(j))<=min(F(grid))+1e-7);
 if xprox(j)~=0,assert(abs(xprox(j)-v(j)+.5*a*sign(xprox(j))/sqrt(abs(xprox(j))))<1e-11);end
end
pass('Global half-power proximal minimizers');
mu=.1;[~,rg]=np_ssqp_smoothing(P,[.3;-.4],mu);fd=zeros(2,1);
for j=1:2,e=zeros(2,1);e(j)=h;fd(j)=(np_ssqp_smoothing(P,[.3;-.4]+e,mu)-np_ssqp_smoothing(P,[.3;-.4]-e,mu))/(2*h);end
assert(norm(fd-rg)<1e-7);pass('SSQP smoothing gradient');
np_work_budget('start',Inf);[y,solved,qmeta]=np_qp(eye(2),[-2;-3],[],[],[1 1],1,zeros(2,1),ones(2,1),[.5;.5],o);
assert(solved,'CRN:QPCheck', ...
 'Native QP failed. Backend: %s\nExit flag: %g\n%s\n%s', ...
 qmeta.backend,qmeta.flag,qmeta.identifier,qmeta.message);
assert(norm(y-[0;1])<1e-7,'CRN:QPCheck', ...
 'QP returned %s; expected [0;1]. Error: %.3e. Backend: %s; exit flag: %g.', ...
 mat2str(y,16),norm(y-[0;1]),qmeta.backend,qmeta.flag);
pass('Constrained QP and bounds');
D=np_generate_regression(cfg.seed,20,[30 20],50,10,.5,.1);
assert(nnz(D.z1)==10&&nnz(D.z2)==10&&nnz(D.z1.*D.z2)==5);
assert(rank(D.C)==5&&size(D.C,1)==6&&norm(D.C*D.truth)<1e-10);pass('Paired regression and redundant affine representation');
Ds=np_generate_svm(cfg.seed,20,40,20,50,10,.3,.05);
assert(nnz(Ds.truth(1:end-1))==10&&all(Ds.truth(Ds.support_union)==1/sqrt(10)));
assert(rank(Ds.C)==5&&norm(Ds.C*Ds.truth)<1e-10);pass('Synthetic classification model');
for split=0:29
 W=np_wdbc(root,cfg.seed+split,cfg);I=[W.train_index;W.validation_index;W.test_index];
 assert(numel(unique(I))==569&&nnz(W.y==1)==127&&nnz(W.yval==1)==42&&nnz(W.ytest==1)==43);
 assert(max(abs(mean(W.X,1)))<1e-12&&max(abs(std(W.X,0,1)-1))<1e-12);
 assert(norm(W.Xval-bsxfun(@rdivide,bsxfun(@minus,W.raw_X(W.validation_index,:),W.training_mean),W.training_sample_std),'fro')<1e-12);
end
pass('All 30 WDBC splits, labels, and training-only standardization');
W1=np_wdbc(root,cfg.seed,cfg);W2=np_wdbc(root,cfg.seed,cfg);assert(isequal(W1,W2));pass('WDBC reproducibility');
for method={'hybrid','no_newton'}
 O=np_solve(P,zeros(2,1),method{1},o);assert(O.success&&O.iterations==0&&O.initial_stationary);
 O=np_solve(P,[1.5;-1.5],method{1},o);audit(O,o);assert(all(isfinite(O.x)));pass(['CRN core: ' method{1}]);
end
P=np_problem('regression',D,[.495 .495 .01],'free');
for method={'CRN','CRN0','IReNA','HpgSRN'}
 switch method{1}
  case 'CRN',O=np_solve(P,D.x0,'hybrid',o);
  case 'CRN0',O=np_solve(P,D.x0,'no_newton',o);
  case 'IReNA',O=np_irena(P,D.x0,o);
  otherwise,O=np_hpgsrn(P,D.x0,o);
 end
 assert(all(isfinite(O.x))&&O.feasibility<=o.feas_tol&&O.counts.total_qp<=o.max_qp);
 if any(strcmp(method{1},{'CRN','CRN0'})),audit(O,o);end
 pass(['Regression execution: ' method{1}]);
end
Ps=np_problem('svm',Ds,[.01 .495 .495],'free');
for penalty={'l1','l2'}
 O=np_convex_svm(Ps,Ds.x0,penalty{1},o);assert(all(isfinite(O.x))&&O.feasibility<=o.feas_tol);
 assert(O.success&&O.residual<=o.tol);pass(['Convex SVM certificate: ' penalty{1}]);
end
O=np_ssqp(Ps,Ds.x0,o);assert(all(isfinite(O.x))&&O.feasibility<=o.feas_tol);
assert(O.counts.total_qp==O.counts.candidate_qp+O.counts.certificate_qp);pass('SSQP execution and QP accounting');
% The convex l1 zero subgradient is bounded; the fractional zero subgradient is free.
Pzero=Ps;xzero=zeros(Ps.d,1);[~,rf]=np_certificate(Pzero,xzero,[],o);
[~,rc]=np_convex_svm_certificate(Pzero,xzero,'l1',o);assert(isfinite(rf)&&isfinite(rc));pass('Separate convex and fractional certificates');
folder=fullfile(root,'results_checks');if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'checks.mat'),'checks','-v7');
f=fopen(fullfile(folder,'checks.txt'),'wt');for j=1:numel(checks),fprintf(f,'PASS: %s\n',checks{j});end;fclose(f);
fprintf('All %d checks passed.\n',numel(checks));
function pass(name)
 checks{end+1}=name;fprintf('PASS: %s\n',name);
end
end
function audit(O,o)
assert(O.feasibility<=o.feas_tol&&O.counts.total_qp<=o.max_qp);
assert(O.counts.total_qp==O.counts.candidate_qp+O.counts.fallback_qp+O.counts.certificate_qp);
for k=1:size(O.history,1)
 if strcmp(O.branches{k},'fallback')
  assert(O.history(k,5)>=o.eta/2-1e-8&&O.history(k,7)<=o.rho+1e-8);
 else
  assert(O.history(k,5)>=o.aN-1e-8&&O.history(k,6)<=o.BN+1e-8);
 end
end
end
