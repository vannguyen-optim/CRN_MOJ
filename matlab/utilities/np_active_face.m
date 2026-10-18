function [B,signature]=np_active_face(P,x)
% Literal original activity only: tolerances never create normal directions.
ia=find(P.A*x==P.b);lo=find(x==P.lb & isfinite(P.lb));hi=find(x==P.ub & isfinite(P.ub));
E=eye(P.d);B=[P.C;P.A(ia,:);E(lo,:);E(hi,:)];
signature=[ia;size(P.A,1)+lo;size(P.A,1)+P.d+hi];
end
