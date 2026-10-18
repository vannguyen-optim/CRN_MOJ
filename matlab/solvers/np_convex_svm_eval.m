function [value,g,regularizer]=np_convex_svm_eval(P,x,penalty)
% ORIGINAL convex predictive model; bias is not penalized.
% Class coefficients exactly match fractional scalarization divided by lambda1.
assert(strcmp(P.kind,'svm')&&P.lambda(1)>0,'Requires a positive-weight SVM problem');
x=x(:);c=P.c/P.lambda(1);r=max(0,1-P.y.*(P.Xa*x));
g=-2*P.Xa'*(c.*P.y.*r);
switch penalty
 case 'l1',regularizer=sum(abs(x(1:end-1)));
 case 'l2',regularizer=.5*sum(x(1:end-1).^2);g(1:end-1)=g(1:end-1)+x(1:end-1);
 otherwise,error('Unknown convex SVM penalty: %s',penalty);
end
value=regularizer+sum(c.*r.^2);
end
