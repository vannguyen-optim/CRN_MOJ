function cfg=CONFIG_EXPERIMENT3()
% Experiment settings inherit the common publication protocol.
cfg=CRN_CONFIG_PUBLICATION();
cfg.scopes={'svm','wdbc'};
end
