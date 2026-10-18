function out=np_l1_svm(P,x,o)
% Convex l1 predictive baseline with matched class coefficients.
out=np_convex_svm(P,x,'l1',o);
end
