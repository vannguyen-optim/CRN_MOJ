function folder=CRN_RUN_SMOKE()
% Small execution check, with its own result directory.
root=CRN_SETUP();cfg=CRN_CONFIG_PUBLICATION();cfg.output=fullfile(root,'results_smoke');
cfg.repetitions=1;cfg.low_grid=3;cfg.reg_dimensions=20;cfg.reg_overlap=.5;
cfg.reg_n=[30 20];cfg.reg_ntest=50;cfg.reg_gamma=.01;
cfg.svm_dimensions=20;cfg.svm_prevalence=.3;cfg.svm_n=[40 20 50];
cfg.zeta=.01;cfg.q=.5;cfg.c=1;cfg.options.max_outer=60;
cfg.options.max_inner=150;cfg.options.max_qp=200;cfg.local_diagnostics=false;
folder=CRN_RUN_PUBLICATION(cfg);
end
