function status=STATUS_PUBLICATION(folder)
root=CRN_SETUP();
if nargin==0,folder=fullfile(root,'results_publication');end
folders=np_result_folders(folder);
status=struct('Experiment',{},'Attempted',{},'Completed',{},'Successful',{},'Pending',{});
for k=1:numel(folders)
 a=load(fullfile(folders{k},'manifest.mat'),'jobs');done=0;success=0;
 for j=1:numel(a.jobs)
  path=fullfile(folders{k},'runs',[a.jobs(j).id '.mat']);
  if exist(path,'file')
   r=load(path,'record');done=done+1;success=success+logical(r.record.out.success);
  end
 end
 [~,name]=fileparts(folders{k});
 status(end+1)=struct('Experiment',name,'Attempted',numel(a.jobs),'Completed',done, ...
  'Successful',success,'Pending',numel(a.jobs)-done);
 fprintf('%s: %d/%d completed; %d successful; %d pending.\n', ...
  name,done,numel(a.jobs),success,numel(a.jobs)-done);
end
end
