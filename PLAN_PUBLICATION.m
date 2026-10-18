function plan=PLAN_PUBLICATION(cfg)
% Count all attempted training models without generating data or solving QPs.
CRN_SETUP();
if nargin==0,cfg=CRN_CONFIG_PUBLICATION();end
np_validate_config(cfg);
plan=struct('Experiment',{},'Models',{});
weights=(cfg.low_grid-1)*(cfg.low_grid-2)/2;
for j=1:numel(cfg.scopes)
 scope=cfg.scopes{j};
 switch scope
  case 'lowdim'
   count=3*size(cfg.low_starts,2)*weights;
  case 'regression'
   methods=3*numel(cfg.reg_geometries)+sum(strcmp(cfg.reg_geometries,'free'));
   count=cfg.repetitions*numel(cfg.reg_dimensions)*numel(cfg.reg_overlap)*numel(cfg.reg_gamma)*methods;
  otherwise
   per_split=(3+2*logical(cfg.fixed_convex_comparisons))*numel(cfg.zeta)*numel(cfg.q) ...
    +2*numel(cfg.c)*numel(cfg.q);
   if strcmp(scope,'svm')
    groups=numel(cfg.svm_dimensions)*numel(cfg.svm_prevalence)*numel(cfg.svm_geometries);
   else
    groups=1;
   end
   count=cfg.repetitions*groups*per_split;
 end
 plan(end+1)=struct('Experiment',np_experiment_folder(scope),'Models',count);
 fprintf('%s: %d models.\n',plan(end).Experiment,count);
end
fprintf('Total: %d attempted models.\n',sum([plan.Models]));
end
