function [t,x]=noisegen(n,samplef)
%this function simulates the actual sampling of noise. The noise 'source'
%is the same, only the sampling can be adjusted

m=1e4;
samplef2=1e4;
t2=linspace(0,(m-1)*1/samplef2,m);
x2=randn(1,m);


z(1:m)=0;
a=0.99;
for i=1:(m-1)
    z(i+1)=a*z(i)+randn(1);
end
x2=z;
%x2=ones(size(x2));

%x2=1*sin(2*pi*1e3*t2);

t=linspace(0,(n-1)*1/samplef,n);
x=resample(x2,samplef,samplef2); %works also as a antialiasing filter

if n>length(x)
    n=length(x);
    disp('too long!')
end

x=x(1:n);
t=t(1:n);