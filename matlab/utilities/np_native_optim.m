function info=np_native_optim()
% Select the MathWorks solvers used by this experiment's recorded protocol.
% MOSEK provides homonymous quadprog/linprog functions with a different API.
q0=which('quadprog');l0=which('linprog');
native_folder=fullfile(matlabroot,'toolbox','optim','optim');
if ~exist(native_folder,'dir')
 error('CRN:OptimizationToolbox', ...
  'MathWorks Optimization Toolbox is required. MATLAB root: %s',matlabroot);
end
if ~is_native(q0)||~is_native(l0)||~is_mathworks_options(which('optimoptions'))
 options_folder=fullfile(matlabroot,'toolbox','shared','optimlib');
 if exist(options_folder,'dir'),addpath(options_folder,'-begin');end
 addpath(native_folder,'-begin');
 rehash;
end
q=which('quadprog');l=which('linprog');options=which('optimoptions');
if ~is_native(q)||~is_native(l)||~is_mathworks_options(options)
 error('CRN:SolverPath', ...
  ['MathWorks solvers are not first on the MATLAB path.\n' ...
   'quadprog: %s\nlinprog: %s\noptimoptions: %s'],q,l,options);
end
if ~license('test','Optimization_Toolbox')
 error('CRN:OptimizationLicense','A MathWorks Optimization Toolbox license is required.');
end
info=struct('quadprog',q,'linprog',l,'optimoptions',options, ...
 'quadprog_before_setup',q0,'linprog_before_setup',l0, ...
 'MATLAB',version,'release',version('-release'));
end
function yes=is_native(file)
prefix=[matlabroot filesep 'toolbox' filesep 'optim' filesep];
if ispc,yes=strncmpi(file,prefix,numel(prefix));
else,yes=strncmp(file,prefix,numel(prefix));end
end

function yes=is_mathworks_options(file)
prefix=[matlabroot filesep 'toolbox' filesep];
if ispc,yes=strncmpi(file,prefix,numel(prefix));
else,yes=strncmp(file,prefix,numel(prefix));end
end
