function [data, header] = read_data(file)
% Read data in TDT SEV file
% 
% INPUT
%    file : string, filename of SEV file
% 
% OUTPUT
%    data : 1D data in same format as in SEV
%    header: struct with header information

% read header
HEADERSIZE = 40; % bytes
header = read_header(file);

% open file
fid = fopen(file, 'rb');
if fid < 0
    error(['Could not open ' file])       
end

% read data
fseek(fid, 40, 'bof');
data = fread(fid, inf, ['*' header.dForm]);

fclose(fid);