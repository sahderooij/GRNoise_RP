function WriteSRONcsv(file,data,header,precision)
% This function writes a file csv file containing the specified data with
% the requested precision. Before the first dataline a line of column
% headers is written (also comma separated) if desired.
%
%INPUT:
%   file == A string that is the complete path of the file. Ending in '\filename.csv'
%   data == An MxN array containing the data that needs to be written to file.
%   header == There are three options for this input parameter:
%             1) Empty [] == no header
%             2) A single string. This will be written as header
%                line.
%             3) A Nx1 cell vector containing strings. These will be
%             concatonated using comma's to create the header string.
%   precision == A string describing the precision used to write the data
%                into the csv. See Matlab help of dlmwrite for details:
%                http://www.mathworks.nl/help/techdoc/ref/dlmwrite.html
%                Suggested Default: '%.e' [scientific, with  digits behind .]
%
% Note: when using header option 3 a warning will be given is the number of
% vector elements is not equal to the number of columns in data.
%
%OUTPUT:
%
%SUBROUTINES:
%
%REQUIRED FILES:
%
%VERSION: 1.0
%
%DATE: August 8, 2012
%AUTHOR: Reinier Janssen
%==========================================================================

%==========================================================================
% CONSTRUCT THE HEADER
%==========================================================================

%Check class of the header input variable.
VariableInfo = whos('header'); %Get properties of the 'header' variable
if strcmp(VariableInfo.class,'cell') %If cell vector
    %Give a warning if there is a discrepancy between header and data size.
    if VariableInfo.size(1) ~= size(data,2)
        fprintf(['Warning (WriteSRONcsv): ',...
'The number of header elements does not correspond to the number of data columns.\n'])
        fprintf('HINT: Check if the variables "data" and "header" have the right orientation\n')
    end
    
    %Concatonate the cells into a single string.
    headerstring = header{1,1};
	for p=2:VariableInfo.size(1)
        headerstring = [headerstring,',',header{p,1}];
    end
else
    %If not a cell array, but just a string (or empty).
    headerstring = header;
end

%==========================================================================
% WRITE THE HEADER
%==========================================================================

%Open the file
fid = fopen(file,'w'); %OVERWRITES EXISTING FILES
if isempty(headerstring) == 0
    %Write the header string
    fprintf(fid,[headerstring,'\n']);
end
%Close the file
fclose(fid);

%==========================================================================
% WRITE THE DATA
%==========================================================================

%write the data with comma separation to the file.
dlmwrite(file,data, '-append','newline', 'pc', 'precision', precision,'delimiter', ',');