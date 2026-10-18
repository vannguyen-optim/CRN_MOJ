function S=np_report(folder)
load(fullfile(folder,'manifest.mat'),'jobs');load(fullfile(folder,'protocol.mat'),'cfg');
rows=struct([]);
for j=1:numel(jobs)
 file=fullfile(folder,'runs',[jobs(j).id '.mat']);
 completed=exist(file,'file')~=0;
 if completed,a=load(file,'record');record=a.record;
 else
  record=struct('job',jobs(j),'out',struct('success',false,'status','pending'));
 end
 J=record.job;O=record.out;
 r=struct('ID',J.id,'Experiment',J.experiment,'Method',J.method,'Geometry',J.geometry, ...
  'Realization',J.realization,'Dimension',J.dimension,'Overlap',J.overlap, ...
  'Prevalence',J.prevalence,'Start',J.start,'Lambda1',J.lambda(1), ...
  'Lambda2',J.lambda(2),'Lambda3',J.lambda(3),'Zeta',J.zeta,'Q',J.q,'C',J.c, ...
  'Role',job_role(J),'Completed',completed,'Success',O.success,'Status',O.status,'Seconds',value(O,'seconds'), ...
  'CPUSeconds',value(O,'cpu_seconds'),'Iterations',value(O,'iterations'),'Residual',value(O,'residual'), ...
  'Feasibility',value(O,'feasibility'),'Objective',value(O,'objective'), ...
  'Selected',value(O,'selected'),'NQP',NaN,'CandidateQP',NaN,'FallbackQP',NaN, ...
  'CertificateQP',NaN,'Nlin',NaN,'ReferenceGap',value(O,'reference_gap'), ...
  'MSE1',value(O,'mse1'),'MSE2',value(O,'mse2'), ...
  'SupportChanges',value(O,'support_changes'),'LastSupportChange',value(O,'last_support_change'), ...
  'ReducedMineig',value(O,'reduced_mineig'),'FaceMineig',value(O,'face_mineig'), ...
  'MinMargin',value(O,'min_margin'),'NativeResidual',value(O,'native_residual'), ...
  'ValidationBA',NaN,'ValidationLoss',NaN,'FixedKey',fixed_key(J),'PredictionKey',prediction_key(J));
 if isfield(O,'counts')
  K=O.counts;r.NQP=value(K,'total_qp');r.CandidateQP=value(K,'candidate_qp');
  r.FallbackQP=value(K,'fallback_qp');r.CertificateQP=value(K,'certificate_qp');r.Nlin=value(K,'linear_solves');
 end
 if isfield(O,'validation'),r.ValidationBA=O.validation.balanced;r.ValidationLoss=O.validation.loss_plus+O.validation.loss_minus;end
 rows=[rows;r];
end
T=text_table(rows);writetable(T,fullfile(folder,'tables','all_runs.csv'));
% Fixed-configuration summaries: never pool dimensions/geometries/weights.
opt=struct([]);keys=unique(T.FixedKey,'stable');
for k=1:numel(keys)
 I=strcmp(T.FixedKey,keys{k});methods=unique(T.Method(I),'stable');
 for j=1:numel(methods)
  if any(strcmp(methods{j},{'l1_SVM','l2_SVM'})),continue;end
  A=T(I&strcmp(T.Method,methods{j}),:);B=A(A.Success,:);
  r=struct('Key',keys{k},'Experiment',A.Experiment{1},'Method',methods{j}, ...
   'Geometry',A.Geometry{1},'Dimension',A.Dimension(1),'Overlap',A.Overlap(1), ...
   'Prevalence',A.Prevalence(1),'Start',A.Start(1),'Lambda1',A.Lambda1(1), ...
   'Lambda2',A.Lambda2(1),'Lambda3',A.Lambda3(1),'Attempted',height(A),'Success',height(B));
  for field={'Iterations','Seconds','Objective','Selected','NQP','Nlin','MSE1','MSE2','ReferenceGap'}
   f=field{1};[av,sd]=stats(B.(f));r.([f 'Mean'])=av;r.([f 'SD'])=sd;
  end
  r.MaximumGap=max_or_nan(B.ReferenceGap);opt=[opt;r];
 end
end
OPT=text_table(opt);writetable(OPT,fullfile(folder,'tables','optimization_by_configuration.csv'));
% Low-dimensional table pools only starts within the explicitly specified diagnostic.
LOW=struct([]);
for diagnostic=0:1
 methods={'CRN','CRN0','SSQP'};
 for j=1:numel(methods)
  I=strcmp(T.Experiment,'lowdim')&strcmp(T.Method,methods{j}) & ((T.Start==3)==logical(diagnostic));
  A=T(I,:);if isempty(A),continue;end;B=A(A.Success,:);
  r=struct('Diagnostic',diagnostic,'Method',methods{j},'Attempted',height(A),'Success',height(B));
  for field={'Iterations','Seconds','ReferenceGap'},f=field{1};[av,sd]=stats(B.(f));r.([f 'Mean'])=av;r.([f 'SD'])=sd;end
  r.MaximumGap=max_or_nan(B.ReferenceGap);LOW=[LOW;r];
 end
end
if ~isempty(LOW),LOW=text_table(LOW);writetable(LOW,fullfile(folder,'tables','experiment1.csv'));end
% Selection uses validation only. Test evaluation begins after an index is chosen.
selected=struct([]);predKeys=unique(T.PredictionKey,'stable');
for k=1:numel(predKeys)
 if strcmp(predKeys{k},'none'),continue;end
 I=strcmp(T.PredictionKey,predKeys{k});methods=unique(T.Method(I),'stable');
 for j=1:numel(methods)
  for split=1:cfg.repetitions
    eligible=find(I&strcmp(T.Method,methods{j})&T.Realization==split&T.Success& ...
     ~strcmp(T.Role,'fixed')&isfinite(T.ValidationBA)&isfinite(T.ValidationLoss)&isfinite(T.Selected));
   r=struct('Key',predKeys{k},'Method',methods{j},'Split',split,'Valid',false, ...
    'ID','','Zeta',NaN,'Q',NaN,'C',NaN,'ValidationBA',NaN,'Seconds',NaN,'Selected',NaN, ...
    'RecallPlus',NaN,'RecallMinus',NaN,'BalancedAccuracy',NaN,'Accuracy',NaN,'LossPlus',NaN,'LossMinus',NaN);
   if ~isempty(eligible)
    parameter=T.Zeta(eligible);if contains(methods{j},'_SVM'),parameter=T.C(eligible);end
    ranking=[-T.ValidationBA(eligible),T.Selected(eligible),T.ValidationLoss(eligible),parameter,T.Q(eligible)];
    [~,order]=sortrows(ranking,[1 2 3 4 5]);index=eligible(order(1));
    chosen=load(fullfile(folder,'runs',[jobs(index).id '.mat']),'record');O=chosen.record.out;
    load(fullfile(folder,jobs(index).input),'D');test=np_class_metrics(D.Xtest,D.ytest,O.x);
    r.Valid=true;r.ID=T.ID{index};r.Zeta=T.Zeta(index);r.Q=T.Q(index);r.C=T.C(index);
    r.ValidationBA=100*T.ValidationBA(index);r.Seconds=O.seconds;r.Selected=O.selected;
    r.RecallPlus=100*test.recall_plus;r.RecallMinus=100*test.recall_minus;
    r.BalancedAccuracy=100*test.balanced;r.Accuracy=100*test.accuracy;
    r.LossPlus=test.loss_plus;r.LossMinus=test.loss_minus;
   end
   selected=[selected;r];
  end
 end
end
SELECT=[];PRED=[];
if ~isempty(selected)
 SELECT=text_table(selected);writetable(SELECT,fullfile(folder,'tables','validation_selected_models.csv'));
 pred=struct([]);keys=unique(SELECT.Key,'stable');
 for k=1:numel(keys)
  methods=unique(SELECT.Method(strcmp(SELECT.Key,keys{k})),'stable');
  for j=1:numel(methods)
   A=SELECT(strcmp(SELECT.Key,keys{k})&strcmp(SELECT.Method,methods{j}),:);B=A(A.Valid,:);
   M=T(find(strcmp(T.PredictionKey,keys{k}),1),:);
   r=struct('Key',keys{k},'Experiment',M.Experiment{1},'Dimension',M.Dimension, ...
    'Prevalence',M.Prevalence,'Geometry',M.Geometry{1},'Method',methods{j},'Attempted',height(A),'Valid',height(B));
   for field={'RecallPlus','RecallMinus','BalancedAccuracy','Accuracy','Selected','LossPlus','LossMinus'}
    f=field{1};[av,sd]=stats(B.(f));r.([f 'Mean'])=av;r.([f 'SD'])=sd;
   end
   pred=[pred;r];
  end
 end
 PRED=text_table(pred);writetable(PRED,fullfile(folder,'tables','test_performance.csv'));
end
% Tuning includes every eligible-parameter training attempt, including failures.
tuning=struct([]);
for k=1:numel(predKeys)
 if strcmp(predKeys{k},'none'),continue;end
 I=strcmp(T.PredictionKey,predKeys{k})&~strcmp(T.Role,'fixed');
 methods=unique(T.Method(I),'stable');
 for j=1:numel(methods)
  for split=1:cfg.repetitions
   A=T(I&strcmp(T.Method,methods{j})&T.Realization==split,:);
   if isempty(A),continue;end
   finiteTime=A.Completed&isfinite(A.Seconds);
   seconds=NaN;if any(finiteTime),seconds=sum(A.Seconds(finiteTime));end
   tuning=[tuning;struct('Key',predKeys{k},'Method',methods{j},'Split',split, ...
    'Candidates',height(A),'Completed',sum(A.Completed),'Successful',sum(A.Success), ...
    'TimedCandidates',sum(finiteTime),'TuningSeconds',seconds)];
  end
 end
end
TUNING=[];
if ~isempty(tuning)
 TUNING=text_table(tuning);writetable(TUNING,fullfile(folder,'tables','tuning_by_split.csv'));
 grouped=struct([]);keys=unique(TUNING.Key,'stable');
 for k=1:numel(keys)
  methods=unique(TUNING.Method(strcmp(TUNING.Key,keys{k})),'stable');
  for j=1:numel(methods)
   A=TUNING(strcmp(TUNING.Key,keys{k})&strcmp(TUNING.Method,methods{j}),:);
   complete=A.Completed==A.Candidates&A.TimedCandidates==A.Candidates;
   [tm,ts]=stats(A.TuningSeconds(complete));
   grouped=[grouped;struct('Key',keys{k},'Method',methods{j},'Splits',height(A), ...
    'Candidates',sum(A.Candidates),'Completed',sum(A.Completed),'Successful',sum(A.Successful), ...
    'FullyTimedSplits',sum(complete),'TuningSecondsMean',tm,'TuningSecondsSD',ts)];
  end
 end
 writetable(text_table(grouped),fullfile(folder,'tables','tuning_summary.csv'));
end
% Matched convex runs retain their own objectives and stationarity criteria.
MATCHED=T(strcmp(T.Role,'fixed')&ismember(T.Method,{'l1_SVM','l2_SVM'}),:);
if ~isempty(MATCHED)
 writetable(MATCHED,fullfile(folder,'tables','matched_convex_comparisons.csv'));
 grouped=struct([]);keys=unique(MATCHED.FixedKey,'stable');
 for k=1:numel(keys)
  methods=unique(MATCHED.Method(strcmp(MATCHED.FixedKey,keys{k})),'stable');
  for j=1:numel(methods)
   A=MATCHED(strcmp(MATCHED.FixedKey,keys{k})&strcmp(MATCHED.Method,methods{j}),:);B=A(A.Success,:);
   r=struct('Key',keys{k},'Experiment',A.Experiment{1},'Method',methods{j}, ...
    'Geometry',A.Geometry{1},'Dimension',A.Dimension(1),'Prevalence',A.Prevalence(1), ...
    'Lambda1',A.Lambda1(1),'Lambda2',A.Lambda2(1),'Lambda3',A.Lambda3(1), ...
    'Attempted',height(A),'Completed',sum(A.Completed),'Success',height(B));
   for field={'Seconds','CPUSeconds','Objective','Selected','NQP','Residual'}
    z=field{1};name=z;if strcmp(z,'Objective'),name='NativeObjective';end
    [av,sd]=stats(B.(z));r.([name 'Mean'])=av;r.([name 'SD'])=sd;
   end
   grouped=[grouped;r];
  end
 end
 writetable(text_table(grouped),fullfile(folder,'tables','matched_convex_summary.csv'));
end
% Pairing requires identical inputs, starting point, geometry, and scalarization.
pairs=struct([]);
for j=find(strcmp(T.Method,'CRN'))'
 J=jobs(j);base=J.id(1:end-numel(J.method));
 for method={'CRN0','IReNA','HpgSRN','SSQP'}
  index=find(strcmp(T.ID,[char(base) method{1}]),1);
  if isempty(index),continue;end
  jointly=T.Success(j)&&T.Success(index);
  r=struct('Key',T.FixedKey{j},'Method',method{1},'Realization',J.realization,'Start',J.start, ...
   'JointSuccess',jointly,'ObjectiveDifference',NaN,'RuntimeRatio',NaN);
  if jointly
   r.ObjectiveDifference=T.Objective(index)-T.Objective(j);
   if isfinite(T.Seconds(j))&&T.Seconds(j)>0&&isfinite(T.Seconds(index))
    r.RuntimeRatio=T.Seconds(index)/T.Seconds(j);
   end
  end
  pairs=[pairs;r];
 end
end
PAIR=[];PAIRSUMMARY=[];
if ~isempty(pairs)
 PAIR=text_table(pairs);writetable(PAIR,fullfile(folder,'tables','paired_runs.csv'));
 grouped=struct([]);keys=unique(PAIR.Key,'stable');
 for k=1:numel(keys)
  methods=unique(PAIR.Method(strcmp(PAIR.Key,keys{k})),'stable');
  for j=1:numel(methods)
   A=PAIR(strcmp(PAIR.Key,keys{k})&strcmp(PAIR.Method,methods{j}),:);B=A(A.JointSuccess,:);
   [dm,ds]=stats(B.ObjectiveDifference);[rm,rs]=stats(B.RuntimeRatio);
   grouped=[grouped;struct('Key',keys{k},'Method',methods{j},'Attempted',height(A),'JointSuccess',height(B), ...
    'ObjectiveDifferenceMean',dm,'ObjectiveDifferenceSD',ds,'RuntimeRatioMean',rm,'RuntimeRatioSD',rs)];
  end
 end
 PAIRSUMMARY=text_table(grouped);writetable(PAIRSUMMARY,fullfile(folder,'tables','paired_summary.csv'));
end
S=struct('all',T,'optimization',OPT,'lowdim',LOW,'selected',SELECT,'prediction',PRED, ...
 'tuning',TUNING,'matched_convex',MATCHED,'paired',PAIRSUMMARY);
save(fullfile(folder,'summary.mat'),'S','-v7');np_tex_tables(folder,S);
end
function k=fixed_key(J)
k=sprintf('%s_d%d_o%g_p%g_%s_start%d_l%.8g_%.8g_%.8g',J.experiment,J.dimension,J.overlap,J.prevalence,J.geometry,J.start,J.lambda);
end
function k=prediction_key(J)
if any(strcmp(J.experiment,{'svm','wdbc'})),k=sprintf('%s_d%d_p%g_%s',J.experiment,J.dimension,J.prevalence,J.geometry);else,k='none';end
end
function role=job_role(J)
if isfield(J,'role'),role=J.role;
elseif any(strcmp(J.method,{'l1_SVM','l2_SVM'})),role='tuning';
else,role='both';end
end
function v=value(O,f)
v=NaN;if isfield(O,f),v=O.(f);end
end
function [m,s]=stats(v)
v=v(isfinite(v));if isempty(v),m=NaN;s=NaN;else,m=mean(v);s=std(v,0);end
end
function v=max_or_nan(v)
if isempty(v),v=NaN;else,v=max(v);end
end

function T=text_table(S)
T=struct2table(S);fields=T.Properties.VariableNames;
for j=1:numel(fields)
 f=fields{j};if ischar(T.(f)),T.(f)=cellstr(T.(f));end
end
end
