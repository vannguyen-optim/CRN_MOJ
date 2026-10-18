function b = np_power_grad(P,x)
b=zeros(P.d,1); J=P.alpha>0 & x~=0;
b(J)=P.alpha(J).*P.p(J).*sign(x(J)).*abs(x(J)).^(P.p(J)-1);
end
