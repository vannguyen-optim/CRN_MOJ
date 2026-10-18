function ref=np_low_reference(lambda)
% Endpoints, zero, and ALL real stationary roots on each sign half-line.
% Numerical polynomial-root reference, not a formal interval-arithmetic certificate.
values=zeros(2,1);z=zeros(2,1);gap=[];
a=.10*lambda(1)+.15*lambda(2)+.05*lambda(3);
c=2*(lambda(1)+lambda(2))-4*lambda(3);
lin=[2*(-lambda(1)+lambda(2));2*(lambda(1)-lambda(2))];
for j=1:2
 cand=[-2;0;2];
 for s=[-1 1]
  rr=roots([4*lambda(3) 0 0 0 c 0 s*lin(j) a/2]);
  good=abs(imag(rr))<=1e-8 & real(rr)>0 & real(rr)<=sqrt(2)+1e-9;
  cand=[cand;s*real(rr(good)).^2];
 end
 % Avoid interpreter-specific temporary indexing in the MATLAB version below.
 v=(lambda(1)+lambda(2))*cand.^2+lin(j)*cand+lambda(1)+lambda(2) ...
   +lambda(3)*(cand.^2-1).^2+a*sqrt(abs(cand));
 [values(j),i]=min(v);z(j)=cand(i);gap=[gap;cand];
end
ref=struct('root_imag_tolerance',1e-8,'root_endpoint_tolerance',1e-9,'value_tolerance',1e-12,'value',sum(values),'point',z,'candidates',gap,'type','polynomial_roots_double_precision');
end
