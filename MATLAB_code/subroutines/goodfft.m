function [f,spectrum]=goodfft(t,x,window,option)
%performs normalised fft on complex x(t)


n=length(t);
if length(x)~=n
    return
end

tmax=t(n);
samplefrequency=1 /  ((tmax-t(1))/(n-1));
%samplefrequency=1/(t(2)-t(1));

f=linspace(0,1-1/n,n)*samplefrequency;



if window==0
    %rectangular
    windowfunction=ones(1,n);
    %cg=1;
end

if window==1
    %hanning
    windowfunction=hann(n)'; %with hann last and first are zero, coherent gain=0.5
    %cg=0.5;
end

if window==2
    %hamming
    windowfunction=hamming(n)';
    %cg=0.54;
end

if window==3
    %blackman
    windowfunction=blackman(n)';
    %cg=0.42;
end

if window==4
    %bartlett (this is useful for autocorrelation). Estimated ACF windowed
    %by bartlett gives true ACF for windowed x[n]. For sure such a function
    %has all ACF properties
    windowfunction=bartlett(n)';
    %cg=0.5;
end

if window==5
    %flat top
    %windowfunction=window(@flattopwin,n)';
    windowfunction=flattop(n)';
    %cg=0.22
end

cg=sum(windowfunction)/n; %coherent gain
enbw=n * sum(windowfunction.^2) / (sum(windowfunction)^2); %equivalent noise bandwidth
windowfactor=1/cg; %; should be 1/cg. But only in the case when x(t) is comparable at every t, i.e. noisy static data. 
%Not when near t=0 there is a lot of signal and at t larger there is no signal like the covariance. There is no point then to try to increase by dividing with the coherent gain. 

%figure
%plot(t,windowfunction.*x)


timeshiftcorr=-1i*2*pi*f.*t(1); %fft function doesn't take time into account, it assumes for n=1 t=0. So fft properties like x=real,even ->X=real,even are lost

spectrum=exp(timeshiftcorr).*fft(windowfunction.*x)/n *windowfactor;



%if option==2
    %SSB output,i.e. f=0... instaed of -.. 0 .. +..
%    f1=spectrum(1:n/2);
