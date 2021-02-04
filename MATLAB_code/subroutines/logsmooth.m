function [fo,so]=logsmooth(fcart,s,ppd)

%remove negative values
idx=find(fcart>0);
fcart=fcart(idx);
s=s(idx);

%ppd=100; %points per decade


%fcart(1)
%fcart(2)
disp([min(fcart),max(fcart)])
f=log10(fcart); %log f
minimum=floor(min(f));
maximum=ceil(max(f));
decades=maximum-minimum;

fkey=linspace(minimum,maximum,ppd*decades); %proposed log f key array


flen=length(f);
fkeylen=length(fkey);


fkeytel=1;
fotel=1;
exitloop=0;
tel=1;
smoothmin=0; %keep at zero, does not work too well

while exitloop==0
    
    ftomatch=fkey(fkeytel);
    ftomatch2=fkey(fkeytel+1);
    df=(ftomatch2-ftomatch)/2;
    
    
    
    smoothtel=0;
    decimated=0;
    exitloop2=0;
    tel2=tel;
    ftot=0;
    stot=0;
    
    while exitloop2==0
        t=f(tel2);
        
        if t>=(ftomatch-df) && t<=(ftomatch+df)
            %start averaging
            ftot=ftot+f(tel2);
            stot=stot+s(tel2);
            smoothtel=smoothtel+1;
            tel=tel2+1;
            if t>=ftomatch && decimated==0
                %in case smoothtel will be < than smoothmin
                fdec=f(tel2);
                sdec=s(tel2);
                decimated=1;
            end          
                    
        end
        
        tel2=tel2+1;
         
        if tel2>flen exitloop2=1;end
        if t>(ftomatch+df) exitloop2=1;end
    
        if (exitloop2==1) && (smoothtel>0)
            % averaging complete, write this to output
            if (smoothtel>smoothmin)
                fo(fotel)=ftot/smoothtel;
                so(fotel)=stot/smoothtel;
            else
                fo(fotel)=fdec;
                so(fotel)=sdec;
            end
            fotel=fotel+1;
        
        end
    
        if exitloop2==1 break;end

    end
      
    
    fkeytel=fkeytel+1;
    
    if tel>=flen exitloop=1;end
    if fkeytel>=fkeylen exitloop=1;end
    
    if exitloop==1 break;end

    

end


fo=10.^fo;