function [bfit, sigmatime]=fitpulseIQ(I,startfit,stopfit,freq)

Ifit=I(startfit:stopfit);
time=linspace(1,length(Ifit),length(Ifit))/freq;
cguess=0;

size(Ifit)
size(time)

cfit = mean(I(end-200:end));
%cfit=mean(I);

s = fitoptions('Method','NonlinearLeastSquares','Startpoint',[1/10 1e-4]);
ftype = fittype(['a*exp(-(x)/b)'],'options', s);
[result bla]=fit(time',Ifit'-cfit,ftype)

afit=result.a
bfit=result.b
%cfit=result.c
result2=confint(result);%confint extracts the 95% confidence intervals from result
sigmatime=(bfit-result2(1,2))/1.96;%convert 95% interval to standard deviation

fitresult=afit.*exp(-time./bfit)+cfit;

figure;
subplot(1,2,1)
plot(time,Ifit,time,fitresult)
xlabel('time(s)');ylabel('amplitude');
subplot(1,2,2)
semilogy(time,Ifit,time,fitresult)
xlabel('time(s)');ylabel('amplitude');

end