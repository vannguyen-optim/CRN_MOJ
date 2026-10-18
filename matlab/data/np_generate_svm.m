function data=np_generate_svm(seed,d,n,nval,ntest,s,prevalence,flips)
rng(seed,'twister');s=min(s,d);J=randperm(d,s);w=zeros(d,1);
w(J)=1/sqrt(s);N=n+nval+ntest;X=randn(N,d);
% Known Gaussian score distribution; split-specific labels need no fitted test statistic.
bias=-sqrt(2)*erfinv(2*(1-prevalence)-1);scores=X*w+bias;
y=ones(N,1);y(scores<0)=-1;flip=rand(N,1)<flips;y(flip)=-y(flip);
if ~any(y(1:n)==1)||~any(y(1:n)==-1),error('Synthetic training has one class');end
r=min(5,d-1);Cw=zeros(r,d);
while rank(Cw)<r,Cw=randn(r,d);Cw=Cw-(Cw*w)*w';end
C=[Cw zeros(r,1)];
data=struct('X',X(1:n,:),'y',y(1:n), 'Xval',X(n+1:n+nval,:), ...
 'yval',y(n+1:n+nval),'Xtest',X(n+nval+1:end,:), 'ytest',y(n+nval+1:end), ...
 'truth',[w;bias],'support_union',J,'C',[C;C(1,:)],'x0',[.1*randn(d,1);0]);
end
