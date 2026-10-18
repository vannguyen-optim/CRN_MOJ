function jobs=np_prepare(cfg,root)
% Generate paired inputs once. Manifest rows correspond to attempted models.
folder=cfg.output;if ~exist(folder,'dir'),mkdir(folder);end
for name={'inputs','runs','tables','figures'},if ~exist(fullfile(folder,name{1}),'dir'),mkdir(fullfile(folder,name{1}));end,end
protocol=cfg;protocol=rmfield(protocol,{'output','parallel','workers'});
np_validate_config(cfg);
source_text=np_code_snapshot(root);
if any(strcmp(cfg.scopes,'wdbc'))
 protocol.wdbc_source=fileread(fullfile(root,'data','wdbc.data'));
end
meta=fullfile(folder,'protocol.mat');
if exist(meta,'file')
 old=load(meta);assert(isequaln(old.protocol,protocol)&&isequal(old.source_text,source_text), ...
  'This folder belongs to a different configuration/source version. Choose cfg.output anew.');
else
 environment=struct('MATLAB',version,'release',version('-release'),'computer',computer, ...
  'timestamp',datestr(now,30),'toolboxes',{ver},'parallel_requested',cfg.parallel,'workers',cfg.workers);
 save(meta,'cfg','protocol','source_text','environment','-v7');
end
jobs=struct('id',{},'input',{},'experiment',{},'method',{},'lambda',{},'geometry',{}, ...
 'realization',{},'dimension',{},'overlap',{},'prevalence',{},'start',{},'zeta',{},'q',{},'c',{},'role',{});
if any(strcmp(cfg.scopes,'lowdim'))
 for a=1:cfg.low_grid-2
  for b=1:cfg.low_grid-a-1
   lambda=[a b cfg.low_grid-a-b]/cfg.low_grid;
   for s=1:size(cfg.low_starts,2)
    id=sprintf('low_%02d_%02d_start%d',a,b,s);D=struct('x0',cfg.low_starts(:,s));file=cache(id,D);
    for method={'CRN','CRN0','SSQP'},append(id,file,'lowdim',method{1},lambda,'box',0,2,NaN,NaN,s,NaN,NaN,NaN);end
   end
  end
 end
end
if any(strcmp(cfg.scopes,'regression'))
 for d=cfg.reg_dimensions
  for overlap=cfg.reg_overlap
   for r=1:cfg.repetitions
    seed=cfg.seed+100000+10000*find(cfg.reg_dimensions==d)+1000*find(cfg.reg_overlap==overlap)+r-1;
    id=sprintf('reg_d%d_o%d_r%02d',d,round(10*overlap),r);
    file=fullfile('inputs',[id '.mat']);
    if ~exist(fullfile(folder,file),'file')
     D=np_generate_regression(seed,d,cfg.reg_n,cfg.reg_ntest,cfg.sparsity,overlap,cfg.noise);
     D.seed=seed;D.generator='twister';cache(id,D);
    end
    for gi=1:numel(cfg.reg_geometries)
     geometry=cfg.reg_geometries{gi};methods={'CRN','CRN0','SSQP'};
     if strcmp(geometry,'free'),methods={'CRN','CRN0','IReNA','HpgSRN'};end
     for gamma=cfg.reg_gamma
      lambda=[(1-gamma)/2 (1-gamma)/2 gamma];tag=sprintf('%s_%s_g%.6g',id,geometry,gamma);
      for method=methods,append(tag,file,'regression',method{1},lambda,geometry,r,d,overlap,NaN,0,NaN,NaN,NaN);end
     end
    end
   end
  end
 end
end
if any(strcmp(cfg.scopes,'svm'))
 for d=cfg.svm_dimensions
  for pi=cfg.svm_prevalence
   for r=1:cfg.repetitions
    seed=cfg.seed+200000+10000*find(cfg.svm_dimensions==d)+1000*find(cfg.svm_prevalence==pi)+r-1;
    id=sprintf('svm_d%d_p%d_r%02d',d,round(10*pi),r);file=fullfile('inputs',[id '.mat']);
    if ~exist(fullfile(folder,file),'file')
     D=np_generate_svm(seed,d,cfg.svm_n(1),cfg.svm_n(2),cfg.svm_n(3),cfg.sparsity,pi,cfg.flips);
     D.seed=seed;D.generator='twister';cache(id,D);
    end
    for gi=1:numel(cfg.svm_geometries),classify(id,file,'svm',cfg.svm_geometries{gi},r,d,pi);end
   end
  end
 end
end
if any(strcmp(cfg.scopes,'wdbc'))
 for r=1:cfg.repetitions
  id=sprintf('wdbc_r%02d',r);file=fullfile('inputs',[id '.mat']);
  if ~exist(fullfile(folder,file),'file'),D=np_wdbc(root,cfg.seed+r-1,cfg);cache(id,D);end
  classify(id,file,'wdbc','free',r,30,NaN);
 end
end
save(fullfile(folder,'manifest.mat'),'jobs','-v7');
function file=cache(id,D)
 file=fullfile('inputs',[id '.mat']);if ~exist(fullfile(folder,file),'file'),save(fullfile(folder,file),'D','-v7');end
end
function classify(id,file,experiment,geometry,r,d,pi)
 for zeta=cfg.zeta
  for q=cfg.q
   lambda=[zeta (1-zeta)*q (1-zeta)*(1-q)];tag=sprintf('%s_%s_z%.6g_q%.6g',id,geometry,zeta,q);
   for method={'CRN','CRN0','SSQP'},append(tag,file,experiment,method{1},lambda,geometry,r,d,NaN,pi,0,zeta,q,NaN);end
   if cfg.fixed_convex_comparisons
    for method={'l1_SVM','l2_SVM'}
     append([tag '_matched'],file,experiment,method{1},lambda,geometry,r,d,NaN,pi,0,zeta,q,NaN);
     jobs(end).role='fixed';
    end
   end
  end
 end
 for c=cfg.c
  for q=cfg.q
   lambda=[1 c*q c*(1-q)];tag=sprintf('%s_%s_c%.6g_q%.6g',id,geometry,c,q);
   for method={'l1_SVM','l2_SVM'},append(tag,file,experiment,method{1},lambda,geometry,r,d,NaN,pi,0,NaN,q,c);end
  end
 end
end
function append(id,file,experiment,method,lambda,geometry,r,d,overlap,pi,start,zeta,q,c)
 jobs(end+1)=struct('id',[id '_' method],'input',file,'experiment',experiment, ...
  'method',method,'lambda',lambda,'geometry',geometry,'realization',r,'dimension',d, ...
  'overlap',overlap,'prevalence',pi,'start',start,'zeta',zeta,'q',q,'c',c,'role','both');
 if any(strcmp(method,{'l1_SVM','l2_SVM'})),jobs(end).role='tuning';end
end
end
