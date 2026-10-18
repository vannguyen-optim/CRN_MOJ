function np_plot(folder,S)
if nargin<2,load(fullfile(folder,'summary.mat'),'S');end
load(fullfile(folder,'protocol.mat'),'cfg');
fctx=fopen(fullfile(folder,'figures','figure_configurations.txt'),'wt');
fprintf(fctx,'Residual display floor: 1e-14. Original saved values are retained.\n');
fprintf(fctx,'Experiment 1 residual: lambda=(0.3,0.3,0.4), start=(1.5,-1.5).\n');
fprintf(fctx,'Experiment 2: free geometry, overlap=%g, gamma=%g; residual dimension=%d, realization=1.\n', ...
 cfg.reg_overlap(ceil(numel(cfg.reg_overlap)/2)),cfg.reg_gamma(ceil(numel(cfg.reg_gamma)/2)),cfg.reg_dimensions(ceil(numel(cfg.reg_dimensions)/2)));
fprintf(fctx,'Experiment 3: dimension=%d, prevalence=%g, free geometry; fixed-cost zeta=%g, q=%g.\n', ...
 cfg.svm_dimensions(1),cfg.svm_prevalence(1),cfg.zeta(1),cfg.q(1));fclose(fctx);
T=S.all;methods={'CRN','CRN0','SSQP','IReNA','HpgSRN','l1_SVM','l2_SVM'};
colors=lines(numel(methods));set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',11);
% Residuals at the exact low-dimensional scalarization in the manuscript.
I=strcmp(T.Experiment,'lowdim')&T.Start==1&abs(T.Lambda1-.3)<1e-12&abs(T.Lambda2-.3)<1e-12;
if any(I),residual_plot(I,'experiment1_residual');end
I=strcmp(T.Experiment,'lowdim')&T.Start~=3&T.Success;
if any(I)
 W=unique([T.Lambda1(strcmp(T.Experiment,'lowdim')),T.Lambda2(strcmp(T.Experiment,'lowdim')),T.Lambda3(strcmp(T.Experiment,'lowdim'))],'rows','stable');
 f=figure('Visible','off');hold on;
 for j=1:3,A=T(I&strcmp(T.Method,methods{j}),:);if isempty(A),continue;end
  [~,weightindex]=ismember([A.Lambda1 A.Lambda2 A.Lambda3],W,'rows');
  plot(weightindex,max(A.ReferenceGap,0),'.','Color',colors(j,:),'DisplayName',methods{j});
 end
 xlabel('Weight index');ylabel('Scalarized reference gap');legend('Location','best');grid on;savefigures(f,'experiment1_gap');
 % Objective-space nondominated computed points.
 idx=find(I);F=[];labels={};for k=idx',a=load(fullfile(folder,'runs',[T.ID{k} '.mat']));F(:,end+1)=a.record.out.F;labels{end+1}=T.Method{k};end
 keep=true(1,size(F,2));for k=1:size(F,2),keep(k)=~any(all(bsxfun(@le,F,F(:,k)),1)&any(bsxfun(@lt,F,F(:,k)),1));end
 f=figure('Visible','off');scatter3(F(1,keep),F(2,keep),F(3,keep),20,'filled');xlabel('f_1');ylabel('f_2');zlabel('f_3');grid on;savefigures(f,'experiment1_nondominated');
end
if any(strcmp(T.Experiment,'regression'))
 overlap=cfg.reg_overlap(ceil(numel(cfg.reg_overlap)/2));gamma=cfg.reg_gamma(ceil(numel(cfg.reg_gamma)/2));
 I=strcmp(T.Experiment,'regression')&T.Overlap==overlap&strcmp(T.Geometry,'free')&abs(T.Lambda3-gamma)<1e-12;
 f=figure('Visible','off');hold on;
 for j=[1 2 4 5]
  av=nan(size(cfg.reg_dimensions));sd=av;
  for k=1:numel(cfg.reg_dimensions),v=T.Seconds(I&T.Dimension==cfg.reg_dimensions(k)&strcmp(T.Method,methods{j})&T.Success);if ~isempty(v),av(k)=mean(v);sd(k)=std(v,0);end,end
  errorbar(cfg.reg_dimensions,av,sd,'-o','Color',colors(j,:),'DisplayName',methods{j});
 end
 xlabel('Dimension');ylabel('Time (s)');legend('Location','best');grid on;savefigures(f,'experiment2_time');
 d=cfg.reg_dimensions(ceil(numel(cfg.reg_dimensions)/2));residual_plot(I&T.Dimension==d&T.Realization==1,'experiment2_residual');
end
if ~isempty(S.selected)
 A=S.selected;I=startsWith(A.Key,'svm_')&A.Valid;
 if any(I)
  % Fixed dimension, prevalence, geometry; the other configurations remain in CSV.
  keys=unique(A.Key(I),'stable');I=I&strcmp(A.Key,keys{1});
  f=figure('Visible','off');hold on;
  for j=1:numel(methods)
   B=A(I&strcmp(A.Method,methods{j}),:);if isempty(B),continue;end
   errorbar(mean(B.Selected),mean(B.BalancedAccuracy),std(B.BalancedAccuracy,0), ...
    std(B.BalancedAccuracy,0),std(B.Selected,0),std(B.Selected,0),'o', ...
    'Color',colors(j,:),'MarkerFaceColor',colors(j,:),'DisplayName',strrep(methods{j},'_','-'));
  end
  xlabel('Selected features');ylabel('Test balanced accuracy (%)');legend('Location','best');grid on;savefigures(f,'experiment3_tradeoff');
 end
end
I=strcmp(T.Experiment,'svm')&T.Dimension==cfg.svm_dimensions(1)&T.Prevalence==cfg.svm_prevalence(1)& ...
 strcmp(T.Geometry,'free')&abs(T.Zeta-cfg.zeta(1))<1e-12&abs(T.Q-cfg.q(1))<1e-12;
if any(I)
 f=figure('Visible','off');av=nan(1,3);sd=av;
 for j=1:3,v=T.Seconds(I&strcmp(T.Method,methods{j})&T.Success);if ~isempty(v),av(j)=mean(v);sd(j)=std(v,0);end,end
 bar(1:3,av);hold on;errorbar(1:3,av,sd,'.k','LineWidth',1);xticks(1:3);xticklabels(methods(1:3));ylabel('Time (s)');grid on;savefigures(f,'experiment3_time');
end
function residual_plot(I,name)
 if ~any(I),return;end
 f=figure('Visible','off');hold on;
 for j=1:numel(methods)
  idx=find(I&strcmp(T.Method,methods{j}),1);if isempty(idx),continue;end
  file=fullfile(folder,'runs',[T.ID{idx} '.mat']);if ~exist(file,'file'),continue;end
  a=load(file);O=a.record.out;
  if ~isfield(O,'history')||isempty(O.history)||~isfield(O,'elapsed'),continue;end
  y=O.history(:,3);t=O.elapsed(2:end)+O.model_seconds;
  if isfield(O,'initial_residual'),y=[O.initial_residual;y];t=[O.model_seconds;t(:)];end
  semilogy(t,max(y,1e-14),'-','Color',colors(j,:),'DisplayName',methods{j});
  k=O.last_support_change;
  if k>0&&k<=size(O.history,1)
   semilogy(O.elapsed(k+1)+O.model_seconds,max(O.history(k,3),1e-14),'s','Color',colors(j,:),'HandleVisibility','off');
  end
 end
 yline(cfg.options.tol,'--k','HandleVisibility','off');xlabel('Elapsed time (s)');ylabel('Original stationarity residual');legend('Location','best');grid on;savefigures(f,name);
end
function savefigures(f,name)
 exportgraphics(f,fullfile(folder,'figures',[name '.pdf']),'ContentType','vector');
 exportgraphics(f,fullfile(folder,'figures',[name '.png']),'Resolution',300);close(f);
end
end
