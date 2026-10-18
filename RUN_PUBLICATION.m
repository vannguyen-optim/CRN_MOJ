function folder=RUN_PUBLICATION(cfg)
% Compatibility entry point; the CRN runner uses unique internal names.
if nargin==0,folder=CRN_RUN_PUBLICATION();else,folder=CRN_RUN_PUBLICATION(cfg);end
end
