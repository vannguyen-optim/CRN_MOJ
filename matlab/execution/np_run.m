function folder=np_run(cfg,root)
assert(exist('quadprog','file')==2,'Optimization Toolbox with quadprog is required.');
assert(license('test','Optimization_Toolbox'),'Optimization Toolbox license is required.');
folder=cfg.output;jobs=np_prepare(cfg,root);pending=false(numel(jobs),1);
for j=1:numel(jobs),pending(j)=~exist(fullfile(folder,'runs',[jobs(j).id '.mat']),'file');end
ids=find(pending);fprintf('%d attempted models; %d pending.\n',numel(jobs),numel(ids));
if cfg.parallel&&~isempty(ids)
 assert(license('test','Distrib_Computing_Toolbox'),'Parallel Computing Toolbox is required for cfg.parallel=true.');
 pool=gcp('nocreate');if isempty(pool),parpool('Processes',cfg.workers);end
 parfor j=1:numel(ids),np_execute(jobs(ids(j)),cfg);end
else
 for j=1:numel(ids)
  np_execute(jobs(ids(j)),cfg);
  if j==1||mod(j,25)==0||j==numel(ids),fprintf('Completed %d/%d pending models.\n',j,numel(ids));end
 end
end
fprintf('Results: %s\n',folder);
end
