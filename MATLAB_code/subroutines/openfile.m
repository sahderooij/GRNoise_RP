function [filestring]=openfile(filename)
%open a file and read its data line for line: filename in, filestring out

    %get size
    t=dir(filename);
    bytelength=t.bytes;

    [fid,message]=fopen(filename,'r');


    %fscanf(fid,'%s %s')
    exit=0;
    
    %predefine filestring for faster reading of file
    filestring=blanks(bytelength);charspace=1;
    
    while exit==0
    line=fgets(fid);
    if line==-1
        exit=1;
    else
        linelength=length(line);
        %filestring=[filestring line];
        filestring(charspace:charspace+linelength-1)=line;charspace=charspace+linelength;
    end

    end
    fclose(fid);
    
    
    %filestring=filestring(1:charspace-1); %remove the last things, actually not necessary, but for future tinkering
    

