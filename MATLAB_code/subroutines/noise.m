function noise

n=1e3;
samplef=1e3; %Hz because the noise we generate artificially is not related to the number here, the PSD depends on the samplef
%t=linspace(0,(n-1)*1/samplef,n);

aver=10;
for tel=1:aver

%x=0*sin(2*pi*100*t)+randn(1,n);
[t,x]=noisegen(n,samplef);
n=length(x);


%x=z;
standarddeviation=std(x)

    Rxx(tel,:)=xcorr(x,'unbiased'); % unbiased: expectation of ACF, Rxx(0)=1 for white noise. Biased: true ACF of windowed x[n], rectangular window applied at x[n], suppresses the ACF at high lags to smooth the periodogram.
    %Rzz(tel,:)=xcorr(z,'unbiased');   


end

scorr=0;
for tel=1:aver
    %[t1,t2]=correlation(x,x,samplef,10,2);
    if tel==1
        %fcorr=t1;
        %scorr=t2;
        %Rz=Rzz(tel,:);
        Rx=Rxx(tel,:);
        
    else
        %scorr=scorr+t2;
        %Rz=Rz+Rzz(tel,:);
        Rx=Rx+Rxx(tel,:);
    end
end
%Rz=Rz/aver;
Rx=Rx/aver;
%scorr=scorr/aver;

%Rzz=Rz;
Rxx=Rx;

Rvector=(-(n-1):n-1)*1/samplef;

%R now averaged and known

[ff,pp]=realfft(Rvector,Rxx,3);
psd=pp*n/samplef; % times n, because realfft divides by n, divided by samplef because PSD is defined at 1 Hz bandwidth.

disp('Parseval:')
sum(Rxx.^2)
sum(psd.^2)

figure
subplot(2,2,1)

plot(t,x,'o-')
title('x(t)')

subplot(2,2,2)

stem(Rvector,Rxx)
title('autocorrelation')


[f,s]=realfft(t,x,0);
%[f,s]=SSB2DSB(f,s);

subplot(2,2,3)
plot(f,abs(s))
title('fft of x(t)')

subplot(2,2,4)

%[pxx,f]=pwelch(x,hann(length(x)),0,n,samplef,'twosided');
[pxx,f]=pwelch(x,[],[],[],samplef,'onesided');



%the onesided FFT of the autocorrelation
[ff,psd]=SSB2DSB(ff,psd);


[f0,p0]=realfft(Rvector,Rxx,0);
p0=p0*n/samplef;
[f0,p0]=SSB2DSB(f0,p0);


[f3,p3]=realfft(Rvector,Rxx,3);
p3=p3*n/samplef;
[f3,p3]=SSB2DSB(f3,p3);

[f4,p4]=realfft(Rvector,Rxx,4);
p4=p4*n/samplef;
[f4,p4]=SSB2DSB(f4,p4);

[f5,p5]=realfft(Rvector,Rxx,5);
p5=p5*n/samplef;
[f5,p5]=SSB2DSB(f5,p5);

plot(f0,p0,f3,p3,f4,p4,f5,p5) %,f,f.^(-1),f,f.^(-2))
%plot(Rvector2,Rxx2)
%hold off
legend('0','3','4','5')



title('psd of x(t)')
