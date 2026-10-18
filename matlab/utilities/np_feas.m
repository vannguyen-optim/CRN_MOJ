function r = np_feas(P,x)
r=max([0; P.lb-x; x-P.ub; abs(P.C*x-P.rhs); P.A*x-P.b]);
if any(~isfinite(x)),r=Inf;end
end
