function secondspectrum(x,freq,average)
%divide x in average parts, compute SSB PSD of each part, and calculate PSD of
%first PSD

% chop up in average parts
n=length(x);


m=n/average;
floor(m);


%psdsize=fix(m/2)+1;
aver2=100;
cpsdmagmatrix=zeros(average,m/aver2);

for tel=1:average
    x2=x( (tel-1)*m+1:(tel)*m );
    [f1,cpsdmag,cpsdph]=correlation(x2,x2,freq,aver2,2); %no averaging
    cpsdmagmatrix(tel,:)=cpsdmag;    %each row = PSD of each part
end





%cpsdmagmatrix: row contains all f's, second row is same for t+dt, with
%dt=1/ (freq/average)

secondspectrum=zeros(average,length(f1));

for tel=1:length(f1)
    x=cpsdmagmatrix(:,tel)'; %this is S(f) at one f for different times
    [f2,cpsdmag,cpsdph]=correlation(x,x,freq/average,1,2);
    secondspectrum(:,tel)=cpsdmag'; %every column contains now S(f2)
end


%plot all PSDs
figure
subplot(2,1,1)
hold on
for tel=1:average
    plot(f1,10*log10(cpsdmagmatrix(tel,:)));
end
hold off

%remove zero element of f2
secondspectrum=secondspectrum(2:length(f2),:);
f2=f2(2:length(f2));

%remove zero element of f1
secondspectrum=secondspectrum(:,2:length(f1));
f1=f1(2:length(f1));

subplot(2,1,2)
surf(f1,f2,10*log10(secondspectrum)) %    Note that x corresponds to the columns of Z and y corresponds to the rows.
xlabel('f')
ylabel('f2')

return