function calculatepsd

%opens a binary file made by labview and reads 64 bit floating point data.
%labview 'header': First 4 bytes: rows in uint, second 4 bytes: columns in uint. 
%Than data row by row in BIG ENDIAN
%PWELCH = PSD * 2/samplefreq!!!





%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-1k-10m-70-500-on-IQ.bin';
%filename='e:\kids-measurements\TaNb\150nmTa-B4-C60-SiOx2\LEDnoise\k39lednoisenicest\nK39-noled-70-750-on-IQ.bin';

filename='e:\kids-jochem\taled\K9\KID9_64dBm_LED2mA_TDmedT0.099996.bin';

[fid,message]=fopen(filename,'r','ieee-be.l64'); %open in Big endian 64 bit floating point number
message
filenameforsize=dir(filename);



%rows=fread(fid,1,'uint')
%columns=fread(fid,1,'uint')

%fseek(fid,8,-1);



jochem=1;

if jochem==1
    rows=filenameforsize.bytes/8/3;
else
    rows=filenameforsize.bytes/8/2;
end

rows

columns=2
freq=250e3;


%reduces memory
readlength=10000; %read readlength columns at a time
pp=ceil(rows/readlength);
I(1:rows)=0;
Q=I;

idx=1;
for tel=1:pp
    rlen=rows-(tel-1)*readlength;
    if rlen>readlength rlen=readlength;end
    if jochem==1
        Matrix=fread(fid,[columns+1,rlen],'float64')';
        I2=Matrix(:,2)';
        Q2=Matrix(:,3)';
        freq=round(1/abs(Matrix(1,1)-Matrix(2,1)));
    else
        Matrix=fread(fid,[columns,rlen],'float64')';
        I2=Matrix(:,1)';
        Q2=Matrix(:,2)';
    end
    I( idx : idx+rlen-1 )=I2;
    Q( idx : idx+rlen-1 )=Q2;
    idx=idx+rlen;
end

fclose(fid);
clear Matrix;

freq




if jochem==1
    %renormalize IQ
    r=sqrt(I.^2+Q.^2);
    p=atan2(Q,I);
    r=r/mean(r);
    I=r.*cos(p);
    Q=r.*sin(p);
end



%resample to conserve memory
%newfreq=freq/1; %new frequency of resampled data
%I=resample(I,newfreq,freq);
%Q=resample(Q,newfreq,freq);
%rows=length(I);
%freq=newfreq;


t=linspace(0,1,rows)*rows/freq;


figh=figure;
subplot(2,2,1)
plot(t,I,t,Q)
title('I(t) & Q(t)')

r=sqrt(I.^2+Q.^2);
p=atan2(Q,I);
p=mod(p,2*pi);

n=length(I);
p=p-mean(p);
r=r-mean(r);



average=100;

[a,f]=pwelch(I-mean(I),n/average,[],n/average,freq,'onesided');

[b,f]=pwelch(Q-mean(Q),n/average,[],n/average,freq,'onesided');
%[prr,fm]=psd(Q-mean(Q),n/average,freq,n/average);%prr=prr*2/freq;

[c,f2]=csd(I-mean(I),Q-mean(Q),n/average,freq,n/average);  % csd is differently defined than in the book, angle is negative. I use: R_{xy}(tau)=E{ x(t)y(t-tau) }
c=c*2/freq;

if f~=f2
    disp('f and f2 not equal')
end

subplot(2,2,2)

psdII=10*log10(a);
psdQQ=10*log10(b);
psdIQ=10*log10(abs(c));
psdIQphase=-angle(c); % csd is differently defined than in the book, angle is negative. I use: R_{xy}(tau)=E{ x(t)y(t-tau) }

%[ff,pp]=correlation(p,p,freq,5);
figure(figh);

semilogx(f,psdII,f,psdQQ,f2,psdIQ)
title('psd I, Q and csd IQ')

subplot(2,2,3)
semilogx(f2,psdIQphase)
title('phase of cross psd')

writegraph([f';psdII';psdQQ';psdIQ';psdIQphase'],'f II QQ IQMag IQph','psd')