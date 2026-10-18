function np_tex_tables(folder,S)
% Complete tables, one fixed configuration per block; short captions.
if ~isempty(S.lowdim)
 f=open('table_van_oct_27.tex');begin(f,'Stationarity and reference gaps.','van_oct_27','lccccc', ...
  'Method & Success & Iterations & Time (s) & Mean gap & Maximum gap');
 for diagnostic=0:1
  text='Nonzero starting points';if diagnostic,text='Zero-start diagnostics';end
  fprintf(f,'\\multicolumn{6}{c}{%s}\\\\\\hline\n',text);
  A=S.lowdim(S.lowdim.Diagnostic==diagnostic,:);
  for j=1:height(A)
   fprintf(f,'%s & %d/%d & %s & %s & %s & %s\\\\\n',label(A.Method{j}),A.Success(j),A.Attempted(j), ...
    pm(A.IterationsMean(j),A.IterationsSD(j)),pm(A.SecondsMean(j),A.SecondsSD(j)), ...
    pm(A.ReferenceGapMean(j),A.ReferenceGapSD(j)),number(A.MaximumGap(j)));
  end
 end
 finish(f);fclose(f);
end
for experiment={'regression','svm','wdbc'}
 I=strcmp(S.optimization.Experiment,experiment{1});A=S.optimization(I,:);if isempty(A),continue;end
 switch experiment{1}
  case 'regression',name='table_van_oct_28.tex';id='van_oct_28';cap='Sparse regression results.';
   cols='lccccccc';header='Method & Success & Iterations & Time (s) & $\Phi_\lambda$ & Selected & $\mathrm{MSE}_1$ & $\mathrm{MSE}_2$';
  case 'svm',name='table_van_oct_29.tex';id='van_oct_29';cap='Optimization performance for synthetic classification.';
   cols='lcccccc';header='Method & Success & Iterations & Time (s) & $\Phi_\lambda$ & $N_{\rm QP}$ & $N_{\rm lin}$';
  otherwise,name='table_van_oct_35.tex';id='van_oct_35';cap='Optimization performance on WDBC.';
   cols='lcccccc';header='Method & Success & Iterations & Time (s) & $\Phi_\lambda$ & $N_{\rm QP}$ & $N_{\rm lin}$';
 end
 f=open(name);keys=unique(A.Key,'stable');
 for k=1:numel(keys)
  B=A(strcmp(A.Key,keys{k}),:);context(f,B(1,:));begin(f,cap,suffixed(id,k),cols,header);
  for j=1:height(B)
   fprintf(f,'%s & %d/%d & %s & %s & %s',label(B.Method{j}),B.Success(j),B.Attempted(j), ...
    pm(B.IterationsMean(j),B.IterationsSD(j)),pm(B.SecondsMean(j),B.SecondsSD(j)),pm(B.ObjectiveMean(j),B.ObjectiveSD(j)));
   if strcmp(experiment{1},'regression')
    fprintf(f,' & %s & %s & %s\\\\\n',pm(B.SelectedMean(j),B.SelectedSD(j)),pm(B.MSE1Mean(j),B.MSE1SD(j)),pm(B.MSE2Mean(j),B.MSE2SD(j)));
   else,fprintf(f,' & %s & %s\\\\\n',pm(B.NQPMean(j),B.NQPSD(j)),pm(B.NlinMean(j),B.NlinSD(j)));end
  end
  finish(f);
 end
 fclose(f);
end
if ~isempty(S.prediction)
 for experiment={'svm','wdbc'}
  A=S.prediction(startsWith(S.prediction.Key,experiment{1}),:);if isempty(A),continue;end
  if strcmp(experiment{1},'svm'),name='table_van_oct_30.tex';id='van_oct_30';cap='Test performance for synthetic classification.';
  else,name='table_van_oct_34.tex';id='van_oct_34';cap='Test performance on WDBC.';end
  f=open(name);keys=unique(A.Key,'stable');
  for k=1:numel(keys)
   B=A(strcmp(A.Key,keys{k}),:);fprintf(f,'%% Configuration: %s\n',keys{k});
   fprintf(f,'\\noindent $d=%d$, geometry: %s',B.Dimension(1),B.Geometry{1});
   if isfinite(B.Prevalence(1)),fprintf(f,', $\\pi_+=%g$',B.Prevalence(1));end
   fprintf(f,'.\n');
   begin(f,cap,suffixed(id,k),'lcccccccc', ...
    'Method & Valid splits & Recall $+$ & Recall $-$ & Balanced accuracy & Accuracy & Selected & $L_+$ & $L_-$');
   for j=1:height(B)
    fprintf(f,'%s & %d/%d',label(B.Method{j}),B.Valid(j),B.Attempted(j));
    for field={'RecallPlus','RecallMinus','BalancedAccuracy','Accuracy','Selected','LossPlus','LossMinus'}
     z=field{1};fprintf(f,' & %s',pm(B.([z 'Mean'])(j),B.([z 'SD'])(j)));
    end
    fprintf(f,'\\\\\n');
   end
   finish(f);
  end
  fclose(f);
 end
end
function f=open(name)
f=fopen(fullfile(folder,'tables',name),'wt');assert(f>=0);fprintf(f,'%% Generated from MATLAB run records.\n');
end
end
function context(f,A)
fprintf(f,'\n\\noindent $d=%d$, geometry: %s',A.Dimension,A.Geometry{1});
if isfinite(A.Overlap),fprintf(f,', overlap: $%g$',A.Overlap);end
if isfinite(A.Prevalence),fprintf(f,', $\\pi_+=%g$',A.Prevalence);end
fprintf(f,', $\\lambda=(%g,%g,%g)$.\n',A.Lambda1,A.Lambda2,A.Lambda3);
end
function begin(f,caption,id,cols,header)
fprintf(f,'\\begin{table}[htbp]\n\\centering\n\\caption{%s}\n\\label{%s}\n\\footnotesize\n\\renewcommand{\\arraystretch}{1.1}\n\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{%s}\n\\hline\n%s\\\\\\hline\n',caption,id,cols,header);
end
function finish(f)
fprintf(f,'\\hline\n\\end{tabular}\n\\end{table}\n\\FloatBarrier\n\n');
end
function s=pm(a,b)
if ~isfinite(a)||~isfinite(b),s='--';else,s=sprintf('$%.4g\\pm%.3g$',a,b);end
end
function s=number(a)
if ~isfinite(a),s='--';else,s=sprintf('$%.4g$',a);end
end
function s=label(a)
s=a;if strcmp(a,'l1_SVM'),s='$\ell_1$-SVM';elseif strcmp(a,'l2_SVM'),s='$\ell_2$-SVM';end
end
function s=suffixed(id,k)
s=id;if k>1,s=sprintf('%s_cfg%d',id,k);end
end
