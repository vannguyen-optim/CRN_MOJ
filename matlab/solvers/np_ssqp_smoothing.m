function [value,g,B]=np_ssqp_smoothing(P,x,mu)
% SSQP smoothing and diagonal upper curvature, equations (4.1), (5.6).
J=find(P.alpha>0);scale=P.alpha(J).^2;t=[scale.*x(J);-scale.*x(J)];m=numel(J);
theta=.5*mu*ones(size(t));I=t>=0&t<=mu;L=t>mu;
theta(I)=t(I).^2/(2*mu)+mu/2;theta(L)=t(L);
prime=zeros(size(t));prime(I)=.5*t(I)./(mu*sqrt(theta(I)));prime(L)=.5./sqrt(t(L));
kappa=2*mu^(-1.5)*double(t>=-mu&t<=2*mu);
g=zeros(P.d,1);B=g;g(J)=scale.*(prime(1:m)-prime(m+1:end));
B(J)=scale.^2.*(kappa(1:m)+kappa(m+1:end));value=sum(sqrt(theta));
end
