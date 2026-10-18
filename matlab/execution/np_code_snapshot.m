function source_text=np_code_snapshot(root)
% Capture all delivered MATLAB source, including experiment configurations.
dirs=[{root},strsplit(genpath(fullfile(root,'matlab')),pathsep), ...
 strsplit(genpath(fullfile(root,'experiments')),pathsep)];
dirs=unique(dirs(~cellfun('isempty',dirs)),'sorted');source_text=cell(0,2);
for k=1:numel(dirs)
 files=dir(fullfile(dirs{k},'*.m'));
 for j=1:numel(files)
  path=fullfile(dirs{k},files(j).name);
  source_text(end+1,:)={strrep(path,[root filesep],''),fileread(path)};
 end
end
end
