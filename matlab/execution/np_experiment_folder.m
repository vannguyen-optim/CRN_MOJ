function name=np_experiment_folder(scope)
switch scope
 case 'lowdim',name='experiment1_stationarity';
 case 'regression',name='experiment2_regression';
 case 'svm',name='experiment3_synthetic';
 case 'wdbc',name='experiment3_wdbc';
 otherwise,error('CRN:Scope','Unknown experiment scope: %s',scope);
end
end
