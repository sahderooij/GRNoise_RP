function [data,header]=skipheaderandconvert(mainstring,keyarray,numberofcolumns,maxheaderlength,entertot)
%getting rid of the header, trigger on keyarray
%headermaxlength=2000;
len=min(maxheaderlength,length(mainstring));
header=mainstring(1:len); %should be big enough to contain the header
arraytosearch=abs(header);


for i=1:(length(arraytosearch)-length(keyarray))
    t=arraytosearch(i:i+length(keyarray)-1);
    if t==keyarray
        break
    end
end
headerendidx=(i-1+length(keyarray)); %index at end of keyarray


%find more enters
if entertot>0
    enterarray=[13 10];
    entertel=0;
    for idx=headerendidx:(length(arraytosearch)-headerendidx - length(enterarray))
        t=arraytosearch(idx:idx+length(enterarray)-1);
        if t==enterarray
            entertel=entertel+1;
            if entertel>=entertot
                break
            end
        end
    end
    headerendidx=idx-1+length(enterarray);
end

header=mainstring(1:headerendidx);
mainstringwithoutheader=mainstring(headerendidx+1:length(mainstring));


%data contains the clean data
data=sscanf(mainstringwithoutheader,'%e',[numberofcolumns,inf])';