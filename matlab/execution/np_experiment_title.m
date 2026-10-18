function title=np_experiment_title(scope)
switch scope
 case 'lowdim',title='Experiment 1: stationarity and reference minima';
 case 'regression',title='Experiment 2: sparse regression';
 case 'svm',title='Experiment 3: synthetic classwise classification';
 case 'wdbc',title='Experiment 3: WDBC';
 otherwise,error('CRN:Scope','Unknown experiment scope: %s',scope);
end
end
