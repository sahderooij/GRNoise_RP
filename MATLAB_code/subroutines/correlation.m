function [f,cpsdmag,cpsdph]=correlation(x,y,freq,average,cartorlog)
%R. Barends, Delft, May 2010
%computes the (cross) power spectral density using the covariance of x,y, sample frequency freq. Data will be
%choppped in average pieces
%outputs f in Hz and cpsd in log10
%covariance = correlation - meanx*meany. Better to calculate covariance
%since correlation will look quite bad at the edges if mean is nonzero, this while mathematically Rxy=Cxy+meanx*meany. FFT of
%covariance is mathematically equal to FFT of correlation (i.e. PSD) except for zero
%frequency component (meanx*meany)
%input: ROW VECTORS!

n=length(x);
if length(y)~=n
    disp('not equal in length')
    return
end

% chop up in average parts
m=n/average;
floor(m);

Rxy=zeros(average,2*m-1);

for tel=1:average
    x2=x( (tel-1)*m+1:(tel)*m );
    y2=y( (tel-1)*m+1:(tel)*m );
    %Rxy(tel,:)=xcorr(x2,y2,'unbiased'); % unbiased: expectation of ACF, Rxx(0)=1 for white noise. Biased: true ACF of windowed x[n], rectangular window applied at x[n], suppresses the ACF at high lags to smooth the periodogram.
    Rxy(tel,:)=xcov(x2,y23,'unbiased'); % unbiased: expectation of ACF, Rxx(0)=1 for white noise. Biased: true ACF of windowed x[n], rectangular window applied at x[n], suppresses the ACF at high lags to smooth the periodogram. Bartlett window.

end

Rsize=size(Rxy(1,:));
R=zeros(Rsize);
for tel=1:average
    %Rz=Rz+Rzz(tel,:);
    R=R+Rxy(tel,:);
end
R=R/average;



Rtime=(-(m-1):m-1)*1/freq; %timing of correlation function

[ff,pp]=goodfft(Rtime,R,4); 


ff(2);
ff(length(ff));
psd=pp*m/freq; % times n, because realfft divides by n, divided by samplef because PSD is defined at 1 Hz bandwidth.

%psd is COMPLEX in case of cross correlation and TWOSIDED


%psd is doublesided: make single sided


[ff,psd]=SSB2DSB(ff,psd);


f=ff;

if cartorlog==1
    %logarithmic output
    cpsdmag=10*log10(real(psd)); %in autocorrelation case, real+even -> phase=0, magnitude=real part of fft
else
    %cartesian output
    cpsdmag=real(psd);%abs(psd);
end

%cpsdmag=abs(psd); %in autocorrelation case, real+even -> phase=0, magnitude=real part of fft

cpsdph=angle(psd);
%cpsdph=mod(cpsdph,2*pi); %0..2*pi





return %!!!





%n=1000;
time=(0:n-1)/freq;
%[time,x]=noisegen(n,freq);
%y=x;



cpsdre=10*log10(abs(real(psd)));
cpsdim=10*log10(abs(imag(psd)));

%Rflip=flipdim(R,2);
%evenR=(R+Rflip)/2;
%oddR=(R-Rflip)/2;
%R=oddR;

%length(Rtime)




%plot
figure
'fig'
subplot(2,3,1)

plot(time,x,time,y)
title('x(t)')

subplot(2,3,2)

%stem(Rtime,R)
plot(Rtime,R,'o-');%,Rtime,evenR,Rtime,oddR)
title('correlation,even,odd')


%[f,s]=realfft(t,x,0);


subplot(2,3,3)
%plot(f,abs(s))
%plot(ff,cpsdph)
hold on
semilogx(f,real(psd),'k')
semilogx(f,imag(psd),'r')
hold off

subplot(2,3,4)
%plot(f,abs(s))
%plot(ff,cpsdph)
hold on
semilogx(f,cpsdre,'k')
semilogx(f,cpsdim,'r')
hold off
title('real and imag of psd of x(t)')



subplot(2,3,5)
plot(f,cpsdmag)

title('magnitude of psd of x(t)')

subplot(2,3,6)

hold on
plot(f,cpsdph,'b.')

%t=flipdim(cpsdph,2);
%plot(f,(cpsdph-t),'ro')

hold off
title('phase of psd of x(t)')

%[pxx,f]=pwelch(x,hann(length(x)),0,n,freq,'twosided');
%[pxx,f]=pwelch(x,[],[],[],freq,'onesided');

%hold on

%loglog(f,pxx)

writegraph([time;x;x-mean(x)],'t x x-mean','corr1')
writegraph([Rtime;R],'t R','corr2')
writegraph([f;cpsdmag],'f mag','corr3')
%bdgngkb

%[ppp,fff]=periodogram(x,hann(length(x)),'twosided',[],samplef);


%PXX=10*log10(pxx);
%PSDre=10*log10(cpsdre);
%PSDim=10*log10(cpsdim);


%semilogx(ff,PSDre,ff,PSDim) %,f,f.^(-1),f,f.^(-2))


%plot(ff,real(psd),ff,imag(psd))
%plot(Rvector2,Rxx2)
%hold off



