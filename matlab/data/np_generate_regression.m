function data=np_generate_regression(seed,d,n,ntest,s,overlap,noise)
rng(seed,'twister');s=min(s,floor(d/2));k=round(overlap*s);perm=randperm(d);
S1=perm(1:s);S2=[S1(1:k) perm(s+1:2*s-k)];
z1=zeros(d,1);z2=z1;z1(S1)=1+rand(s,1);z2(S2)=1+rand(s,1);
% Shared active coefficients agree when their supports overlap.
if k>0,z2(S1(1:k))=z1(S1(1:k));end
zstar=(z1+z2)/2;
A1=randn(n(1),d);A2=randn(n(2),d);At1=randn(ntest,d);At2=randn(ntest,d);
data=struct('A1',A1,'A2',A2,'b1',A1*z1+noise*randn(n(1),1), ...
 'b2',A2*z2+noise*randn(n(2),1),'At1',At1,'At2',At2, ...
 'bt1',At1*z1+noise*randn(ntest,1),'bt2',At2*z2+noise*randn(ntest,1), ...
 'z1',z1,'z2',z2,'truth',zstar,'support_union',union(S1,S2));
r=min(5,d-1);C=zeros(r,d);
while rank(C)<r,C=randn(r,d);C=C-(C*zstar)*(zstar')/(zstar'*zstar);end
data.C=[C;C(1,:)]; % Deliberately rank-deficient representation.
data.x0=.1*randn(d,1); % Common dense start; separately projected for each constraint.
end
