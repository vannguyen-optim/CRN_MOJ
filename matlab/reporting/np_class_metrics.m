function m=np_class_metrics(X,y,z)
y=y(:);z=z(:);
assert(size(X,1)==numel(y)&&size(X,2)+1==numel(z),'Classifier dimensions do not agree.');
assert(all(y==1|y==-1)&&any(y==1)&&any(y==-1),'Both binary classes are required.');
assert(all(isfinite(z))&&all(isfinite(X(:))),'Classifier inputs must be finite.');
score=X*z(1:end-1)+z(end);pred=ones(size(y));pred(score<0)=-1;
r=max(0,1-y.*score);I=y==1;J=y==-1;
m=struct('accuracy',mean(pred==y),'recall_plus',mean(pred(I)==1), ...
 'recall_minus',mean(pred(J)==-1),'balanced',.5*(mean(pred(I)==1)+mean(pred(J)==-1)), ...
 'loss_plus',mean(r(I).^2),'loss_minus',mean(r(J).^2));
end
