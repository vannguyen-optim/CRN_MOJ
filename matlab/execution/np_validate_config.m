function np_validate_config(cfg)
% Validate the supplied protocol before creating or resuming run files.
valid={'lowdim','regression','svm','wdbc'};
assert(iscell(cfg.scopes)&&~isempty(cfg.scopes)&&all(ismember(cfg.scopes,valid)), ...
 'CRN:Config','cfg.scopes must contain supported experiment names.');
assert(numel(unique(cfg.scopes))==numel(cfg.scopes),'CRN:Config','Duplicate scopes.');
assert(ischar(cfg.output)&&~isempty(cfg.output),'CRN:Config','Set a result directory.');
assert(cfg.repetitions>=1&&cfg.repetitions==fix(cfg.repetitions),'CRN:Config','Repetitions must be a positive integer.');
assert(cfg.low_grid>=3&&cfg.low_grid==fix(cfg.low_grid),'CRN:Config','low_grid must be an integer at least 3.');
assert(size(cfg.low_starts,1)==2&&size(cfg.low_starts,2)==3&&all(cfg.low_starts(:,3)==0), ...
 'CRN:Config','Use two nonzero starts and a third literal-zero diagnostic start.');
assert(all(abs(cfg.low_starts(:))<=2),'CRN:Config','Low-dimensional starts must lie in the box.');
assert(all(cfg.reg_gamma>0&cfg.reg_gamma<1)&&all(cfg.zeta>0&cfg.zeta<1), ...
 'CRN:Config','All fractional-power scalarization weights must be strictly positive.');
assert(all(cfg.q>0&cfg.q<1)&&all(cfg.c>0),'CRN:Config','q must be in (0,1) and c positive.');
assert(all(ismember(cfg.reg_geometries,{'free','box','affine'}))&& ...
 all(ismember(cfg.svm_geometries,{'free','affine'})),'CRN:Config','Unsupported geometry.');
assert(all(cfg.reg_overlap>=0&cfg.reg_overlap<=1)&& ...
 all(abs(cfg.sparsity*cfg.reg_overlap-round(cfg.sparsity*cfg.reg_overlap))<1e-12), ...
 'CRN:Config','Support overlap must give an integer intersection.');
assert(all(cfg.reg_dimensions>=2*cfg.sparsity)&&all(cfg.svm_dimensions>=cfg.sparsity+6), ...
 'CRN:Config','Dimensions are too small for the prescribed supports and rank-five constraints.');
assert(all(cfg.svm_prevalence>0&cfg.svm_prevalence<1),'CRN:Config','Prevalence must lie in (0,1).');
assert(cfg.noise>=0&&cfg.flips>=0&&cfg.flips<=1,'CRN:Config','Invalid noise or flip probability.');
assert(all(cfg.wdbc_train+cfg.wdbc_val+cfg.wdbc_test==[212 357]), ...
 'CRN:Config','WDBC class counts must partition all observations.');
o=cfg.options;
assert(o.tol>0&&o.feas_tol>0&&o.qp_tol>0&&o.eta>0&&o.rho>0&&o.chi>0&& ...
 o.aN>0&&o.BN>0&&o.mN>0&&o.beta>0&&o.beta<1,'CRN:Config','Invalid solver parameter.');
limits=[o.max_outer o.max_inner o.max_qp o.qp_max_iter o.newton_trials o.explore_steps];
assert(all(limits>=1&limits==fix(limits)),'CRN:Config','Work limits must be positive integers.');
assert(cfg.workers>=1&&cfg.workers==fix(cfg.workers),'CRN:Config','workers must be a positive integer.');
end
