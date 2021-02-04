% binfilefunction2los(filename,filenamecal,rows,freq,jochem,average,method,plotx,savex,savepath,ppd,rejectpeak,nrsigma,smoothtime,usecalibration)
%R. Barends, Delft, May 2010
%PdV Okt 2017
%clear all

%path = '/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/LT125_VIS/Nov6MUX/TestKID6 402 nm 4000nW 10dB/TD_power/'
%path='/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/LT125_VIS/Noisetest laser off/TD_power/'
%path='/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/LT125_VIS/Dark_100mK-lids 03-10-2017/TD_power/'
%path='/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/LT125_VIS/9KIDs laser off/TD_power/'
%path='/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/LT125_VIS/9KIDs laser on 402 4000nW 10dB/TD_power/'
%path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\TD_Powerlaser1_21dB\'; %LT145
%path = 'C:\Users\Kid\Desktop\laser 120mK\5KIDs laser on 673 160nW 20dB\TD_Power\'; %LT139
%path = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT160_Chip12VISmembrane\Noise_vsT\TD_Power\'
path = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\Hf Goddard Chip 1\Noiselaser402nm_17dB\TD_Power\'
path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165Temp\TD_Power\'
path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165Temp\4KIDs laser on 673 160nW 26dB\TD_Power\'
path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165Temp\4KIDS laser on 1545 160nW 24dB\TD_Power\'
path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165Temp\4KIDs laser on 402 14000nW 18dB\TD_Power\'

path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\4KIDs laser on 1545 50nW 20dB\TD_Power\'
path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165_Chip7VISmembraneIDC\4KIDs laser off\TD_Power\'
path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165_Chip7VISmembraneIDC\4KIDs laser on 250mK 1545 50nW 26dB\TD_Power\'
path = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT165_Chip7VISmembraneIDC\2KIDs laser on 402 14000nW 19dB\TD_Power\'
path = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\Hf Goddard Chip 1\Pulselaser_1dB_4uW\TD_Power\'
%path = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\Hf Goddard Chip 1\Noiselaser_32dB\TD_Power\'
%path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\TDbla\'

%path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\LT165Temp\4KIDs laser on 986 100nW 20dB\TD_Power\'

%path = '\\MARS\kid\KIDonSun\AAA_Scratch temp copies\tempPieter\Noise42mK_Popt-2\TD_Power\'
%path = 'C:\Users\Kid\Desktop\test402nm8uW10dB\TD_Power\'


%path='/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/22Novsubset/10KIDs laser on 402 3000nW 10dB/TD_power/'
%path='/Users/pieterv/SurfdriveDelft/Measurementsfrom2017/LT125_VIS/Noistest trigger 402 3 nW 6us/TD_power/'
%path = '/Volumes/Elements/Nov22lighttightLNA/10KIDs noise/TD_power/'
%path = '/Volumes/Elements/LT132_SiO2patch_30KIDs/TD_Power/'
%filename = 'KID19_90dBm__TDmed_TmK120.bin';
%filename = 'KID18_100dBm__TDmed_TmK120.bin';
%filename = 'KID17_104dBm__TDmed_TmK120.bin';
%filename = 'KID14_79dBm__TDvis20_TmK120.bin';
%filename = 'KID13_86dBm__TDvis20_TmK120.bin';
%filename = 'KID11_102dBm__TDmed_TmK120.bin';
%filename = 'KID10_74dBm__TDvis20_TmK120.bin';
%filename = 'KID8_80dBm__TDvis20_TmK120.bin';
%filename = 'KID7_94dBm__TDmed_TmK120.bin';
%filename = 'KID6_104dBm__TDvis0_TmK120.bin';
%filename = 'KID5_66dBm__TDvis20_TmK120.bin';
%filename = 'KID3_79dBm__TDvis0_TmK120.bin';
filename = 'KID4_101dBm__TDvis2_TmK80.bin';
%filename = 'KID19_94dBm__TDmed_TmK120.bin';


filenamepath = [path filename];
%slow
freq=50e3;              %standard is 50e3;
rows=2e6;               %standard is 2e6; (ie sampling for 40 s at 50 kHz)
% %fast
% freq=1e6;              %standard is 1e6;
% rows=2e5;              %standard is 2e5; (ie sampling for 200 ms at 1 MHz
%vis
freq=1e6;                %1e6 or 2e6
rows=1e6;               %1e6 or 2e6

smoothtime = 2e-5;
smoothit=0;
plotx=1;
noiseav=10; %averages for PSD calculation
ppd = 20;

factorAt = 0.65/0.11;
%vs2 is with subtraction of calibration file (in FFT)

%opens a binary file made by labview and reads 64 bit floating point data.
%labview 'header': First 4 bytes: rows in uint, second 4 bytes: columns in uint. 
%Than data row by row in BIG ENDIAN
%PWELCH = PSD * 2/samplefreq!!!

%%%%%%%%%%%FOR DATA-FILE ON RESONANCE
[fid,message]=fopen(filenamepath,'r','ieee-be.l64'); %open in Big endian 64 bit floating point number
if ~isempty(message)
    fprintf([message '\n']);
end
%filenameforsize=dir(filename);

columns=2;
corroff=1;
crosstime=0; %computes cross-correlation of amp and phase in time

%reduces memory
readlength=10000; %read readlength rows at a time
pp=ceil(rows/readlength);
I(1:rows)=0;
Q=I;

idx=1;
for tel=1:pp
    rlen=rows-(tel-1)*readlength;
    if rlen>readlength; rlen=readlength;end
    Matrix=fread(fid,[columns,rlen],'float64')';
    I2=Matrix(:,1)';
    Q2=Matrix(:,2)';
    I( idx : idx+rlen-1 )=I2;
    Q( idx : idx+rlen-1 )=Q2;
    idx=idx+rlen;
end

fclose(fid);
clear Matrix;

makeIQrp=1; %leave this option 1!!!%here variables I and Q are translated in R and P
if makeIQrp==1
    r=sqrt(I.^2+Q.^2);
    R=r/mean(r); %normalize radius to 1
    p=atan2(Q,I);%This phase is with respect to the axis in the direction (I,Q)=(1,0) in stead of (-1,0) therefore next line
    P=pi-mod(p,2*pi);
    %P=tan(P);
end

figure;
plot(I(1:end),Q(1:end));%hold on;plot(I(1200:1500),Q(1200:1500));
figure
plot(P-mean(P))



t=linspace(0,1,rows)*rows/freq;
tplot=linspace(0,t(end),10);%to plot straight lines based on a few points

n=length(R);

if smoothit==1
    nrpointsfilter=round(freq*smoothtime/2); %ie freq/filterfreq
    evencheck=nrpointsfilter/2-floor(nrpointsfilter/2);
    if evencheck==0
        nrpointsfilter=nrpointsfilter+1;
    end
    Imm=I-mean(I);
    Immsmooth=smooth(Imm,nrpointsfilter)';
    Immsmoothorig=Immsmooth;
    Qmm=Q-mean(Q);
    Qmmsmooth=smooth(Qmm,nrpointsfilter)';
    Qmmsmoothorig=Qmmsmooth;
    stdImm=std(Imm);stdImmsmooth1=std(Immsmooth);
    Rsmooth = smooth(R,nrpointsfilter)';
    Psmooth = smooth(P,nrpointsfilter)';
end
RPadded = (P+(1-R)*factorAt)/2;

    %plot time domain data
if plotx==1
    figh=figure;

    %only plot time domain of normalized amplitude, else fig-file gets too
    %large
    subplot(2,1,1)
    %time domain of I - mean I + sigma levels
    plot(t,R);hold on;
    %plot(t,Rsmooth)
    ylabel('KID amplitude')
    legend('time domain data','smoothed')
    title({path,filename})
    xlabel('time (s)')

    subplot(2,1,2)
    %time domain of I - mean I + sigma levels
    plot(t,P);hold on;
    %plot(t,Psmooth)
    ylabel('KID phase (rad)')
    legend('time domain data','smoothed')
    xlabel('time (s)')

 figure
 semilogy((t)*1000,P)
 ylabel('KID phase (rad)')
 legend('time domain data')
 xlabel('time (s)')


end
    
[srrpartraw,~]=pwelch(R,floor(length(R)/noiseav),[],floor(length(R)/noiseav),freq,'onesided');
[spppartraw,fmb]=pwelch(P,floor(length(P)/noiseav),[],floor(length(P)/noiseav),freq,'onesided');
[~,srrsm]=logsmooth(fmb,srrpartraw,ppd);
[fm2,sppsm]=logsmooth(fmb,spppartraw,ppd);

figure
semilogx(fm2,10*log10(srrsm))
hold on;semilogx(fm2,10*log10(sppsm))
xlabel('Frequency (Hz)')
ylabel('PSD (dBc/Hz)')
legend('amplitude','phase')


% if plotx==1
%     %time domain of smoothed data
%     figure(figh)
% 
%     stdplusfilter=zeros(size(tplot));stdplusfilter(:,:)=nrsigma*stdImmsmooth1;
%     stdminfilter=zeros(size(tplot));stdminfilter(:,:)=-nrsigma*stdImmsmooth1;
% 
%     minstdplus=zeros(size(tplot));minstdplus(:,:)=min5stdImm;
%     minstdmin=zeros(size(tplot));minstdmin(:,:)=-min5stdImm;
%     %time domain of smoothed data
%     subplot(2,2,2)
%     plot(t,Immsmoothorig,tplot,stdplusfilter,'-r',tplot,stdminfilter,'-r',tplot,minstdplus,'-g',tplot,minstdmin,'-g');
%     title('I(t) - mean(I) filtered with quasiparticle lifetime')
%     legend('smoothed time domain data',['I filtered +' num2str(nrsigma) '\sigma'],['I filtered -' num2str(nrsigma) '\sigma'])
% 
%     subplot(2,2,4)
% 
%     plot(t,Immsmooth,tplot,stdplusfilter,'-r',tplot,stdminfilter,'-r',tplot,minstdplus,'-g',tplot,minstdmin,'-g');
%     title(['Rejected parts: ' num2str(nrrejected) '/' num2str(average)])
%     legend('smoothed time domain + rejections',['I filtered +' num2str(nrsigma) '\sigma'],['I filtered -' num2str(nrsigma) '\sigma'])
% 
%     Figfiletd=[filename(1:end-4) '-plot.fig'];%use the filename from the bin-file
%     saveas(figh,Figfiletd,'fig')
%     close(figh) %these fig-files are so large that they need to be closed immediately
% end