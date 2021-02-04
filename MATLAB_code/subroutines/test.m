function test
n=1000;
samplef=1000;
t=linspace(0,(n-1)*1/samplef,n);

freq=10;
x=cos(2*pi*freq*t);
%x=randn(1,n);
[f,a]=realfft(t,x,2);
figure

subplot(2,1,2)
hold on
stem(f,real(a))
stem(f,imag(a),'r')
hold off

subplot(2,1,1)
%plot(t,x)
plot(t,x,t,ifft(a)*n)