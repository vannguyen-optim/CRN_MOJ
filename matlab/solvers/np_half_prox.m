function x=np_half_prox(v,a)
% Global proximal minimizer of 1/2*(x-v)^2+a*sqrt(abs(x)).
% At a tie, choose the zero minimizer.
a=a+zeros(size(v));x=zeros(size(v));I=a==0;x(I)=v(I);
I=a>0 & abs(v)>1.5*a.^(2/3);t=abs(v(I));aa=a(I);
r=(aa/4).*(t/3).^(-1.5);r=min(1,max(-1,r));
x(I)=sign(v(I)).*(2*t/3).*(1+cos(2*pi/3-2*acos(r)/3));
end
