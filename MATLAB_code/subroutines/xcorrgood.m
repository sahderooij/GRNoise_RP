function Rxy=xcorrgood(x,y)
%xcorr of matlab and this function are not correct when x or y has nonzero mean.
%Rxy[m]= Sum x[n] y[n-m]*
%unbiased



n=length(x);

if length(y)~=n
    disp('x and y not equal length')
    return
end



Rxy=zeros(1,2*n-1);
for m=-(n-1):n-1

    

xextra=zeros(1,2*n-1);
yextra=zeros(1,2*n-1);

if m>=0
    %index n is the middle
    xindex=1;
    yindex=m+1;
    xextra(xindex:xindex+n-1)=x;
    yextra(yindex:yindex+n-1)=conj(y);
    s=xextra*yextra';
end


if m<0
    %index n is the middle
    xindex=n;
    yindex=n+m;
    xextra(xindex:xindex+n-1)=x;
    yextra(yindex:yindex+n-1)=conj(y);
    s=xextra*yextra';
end

%overlap=n-|m|
Rxy(n+m)=s/(n-abs(m));

end