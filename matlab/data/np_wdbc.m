function D=np_wdbc(root,seed,cfg)
% Original UCI WDBC file is bundled; no download and no fitted test statistic.
fid=fopen(fullfile(root,'data','wdbc.data'),'rt');assert(fid>=0,'Cannot open WDBC data');
c=onCleanup(@()fclose(fid));a=textscan(fid,['%f%s' repmat('%f',1,30)],'Delimiter',',');
ids=a{1};labels=a{2};X=cell2mat(a(3:end));y=-ones(size(ids));y(strcmp(labels,'M'))=1;
assert(all(isfinite(X(:)))&&isequal(size(X),[569 30])&&nnz(y==1)==212&&nnz(y==-1)==357);
rng(seed,'twister');tr=[];va=[];te=[];remaining=cell(2,1);
for j=1:2
 label=3-2*j;I=find(y==label);I=I(randperm(numel(I)));
 tr=[tr;I(1:cfg.wdbc_train(j))];remaining{j}=I(cfg.wdbc_train(j)+1:end);
end
rng(seed+10000,'twister');
for j=1:2
 rem=remaining{j};rem=rem(randperm(numel(rem)));
 va=[va;rem(1:cfg.wdbc_val(j))];te=[te;rem(cfg.wdbc_val(j)+1:end)];
 assert(numel(rem)-cfg.wdbc_val(j)==cfg.wdbc_test(j));
end
assert(numel(tr)==341&&numel(va)==114&&numel(te)==114);
mu=mean(X(tr,:),1);sd=std(X(tr,:),0,1);sd(sd==0)=1;
Z=bsxfun(@rdivide,bsxfun(@minus,X,mu),sd);
rng(seed+20000,'twister');
D=struct('X',Z(tr,:),'y',y(tr),'Xval',Z(va,:),'yval',y(va), ...
 'Xtest',Z(te,:),'ytest',y(te),'x0',[.1*randn(30,1);0], ...
 'C',zeros(0,31),'ids',ids,'raw_X',X,'raw_y',y,'train_index',tr, ...
 'validation_index',va,'test_index',te,'training_mean',mu,'training_sample_std',sd, ...
 'seed',seed,'validation_test_seed',seed+10000,'start_seed',seed+20000, ...
 'generator','twister','source','https://doi.org/10.24432/C5DW2B');
end
