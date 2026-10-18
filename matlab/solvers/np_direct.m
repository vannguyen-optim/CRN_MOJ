function ok=np_direct(P,x,y,u,o)
d=norm(y-x);
ok=np_feas(P,y)<=o.feas_tol && all(isfinite(u)) && ...
 np_difference(P,x,y)+o.aN*d^2<=0 && norm(u)<=o.BN*d;
end
