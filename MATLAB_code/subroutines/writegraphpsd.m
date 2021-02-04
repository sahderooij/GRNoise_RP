function writegraphpsd(b,header,path,fname)
%R. Barends, Delft, May 2010
%syntax: writegraph( vectors in rows, header, filename)
%writes file in the same directory as this program

%for row vectors:  writegraph([x ; y1 ; y2],'header','filename')

%for matrix: writegraph( [matrix'] ,'header','filename') because it writes
%row vectors in columns;
%global writegraphtel header



%writegraphtel=writegraphtel+1;
a=b';
%writes matrix a (vertically for a row vector b)

filename1=path;

%drivename=thisfunctionfullpath(1); %extract the drive
%filename1=[drivename ':\kids\matlab\'];

filename2='' ; %['-' num2str(writegraphtel)];

filename=[filename1 fname filename2 '.txt'];



[fid,message]=fopen(filename,'w');
if fid==-1
    disp('error opening writegraph file')
    disp(message)
    return
end

    
%write the file
asize=size(a);
rows=asize(1);
columns=asize(2);

if fid~=0
    if isempty(header)==0
        if length(header)>0
            text=['\r\n' header '\r\n \r\n'];
            fprintf(fid,text);
        end
    end
    for r=1:rows
        for c=1:columns
            fprintf(fid,'%e ',a(r,c));
        end
        fprintf(fid,' \r\n');
    end
end

fclose(fid);

fprintf(['Results saved in file: ' fname '\r\n'])

return