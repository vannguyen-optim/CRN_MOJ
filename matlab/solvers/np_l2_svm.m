function out=np_l2_svm(P,x,o)
% Convex l2 predictive baseline with matched class coefficients.
out=np_convex_svm(P,x,'l2',o);
end
