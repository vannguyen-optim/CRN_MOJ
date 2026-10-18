function PLOT_RESULTS(folder)
root=CRN_SETUP();
if nargin==0,folder=fullfile(root,'results_publication');end
folders=np_result_folders(folder);
for j=1:numel(folders)
 np_plot(folders{j});
 fprintf('PDF/PNG figures: %s\n',fullfile(folders{j},'figures'));
end
end
