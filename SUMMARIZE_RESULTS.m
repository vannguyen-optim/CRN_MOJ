function S=SUMMARIZE_RESULTS(folder)
root=CRN_SETUP();
if nargin==0,folder=fullfile(root,'results_publication');end
folders=np_result_folders(folder);summaries=cell(size(folders));
for j=1:numel(folders)
 summaries{j}=np_report(folders{j});
 fprintf('CSV and LaTeX tables: %s\n',fullfile(folders{j},'tables'));
end
if numel(summaries)==1,S=summaries{1};else,S=summaries;end
end
