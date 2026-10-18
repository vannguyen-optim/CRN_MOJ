function folder=RUN_EXPERIMENT2(cfg)
% Run or resume this experiment using the publication settings.
CRN_SETUP();
if nargin==0,cfg=CRN_CONFIG_PUBLICATION();end
cfg.scopes={'regression'};
folder=CRN_RUN_PUBLICATION(cfg);
end
