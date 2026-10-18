function folder=CRN_RUN_PUBLICATION(cfg)
% Run or resume every requested experiment into its own result folder.
root=CRN_SETUP();
if nargin==0,cfg=CRN_CONFIG_PUBLICATION();end
np_validate_config(cfg);
folder=cfg.output;
for j=1:numel(cfg.scopes)
 scope=cfg.scopes{j};local=cfg;local.scopes={scope};
 local.output=fullfile(folder,np_experiment_folder(scope));
 fprintf('\n=== %s ===\n',np_experiment_title(scope));
 np_run(local,root);
 if cfg.export_tables
  np_report(local.output);
  fprintf('CSV and LaTeX tables: %s\n',fullfile(local.output,'tables'));
 end
 if cfg.export_figures
  np_plot(local.output);
  fprintf('PDF/PNG figures: %s\n',fullfile(local.output,'figures'));
 end
end
fprintf('\nFinished requested experiments: %s\n',folder);
end
