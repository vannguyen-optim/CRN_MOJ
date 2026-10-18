function root=CRN_SETUP()
% Add the experiment and shared MATLAB folders to the search path.
root=fileparts(mfilename('fullpath'));
addpath(root,genpath(fullfile(root,'matlab')),genpath(fullfile(root,'experiments')));
np_native_optim();
end
