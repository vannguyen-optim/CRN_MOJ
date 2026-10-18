function folders=np_result_folders(folder)
% Accept one experiment directory or the parent of all experiment results.
if exist(fullfile(folder,'manifest.mat'),'file'),folders={folder};return;end
folders={};
for scope={'lowdim','regression','svm','wdbc'}
 candidate=fullfile(folder,np_experiment_folder(scope{1}));
 if exist(fullfile(candidate,'manifest.mat'),'file'),folders{end+1}=candidate;end
end
assert(~isempty(folders),'CRN:Results','No experiment manifests found in %s.',folder);
end
