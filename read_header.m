function header = read_header(file)
% Read header of TDT SEV file.

ALLOWED_FORMATS = {'single','int32','int16','int8','double','int64'};

% open file
fid = fopen(file, 'rb');

if fid < 0
    error(['Could not open ' file])       
end

% create and fill header struct
header = [];

header.fileSizeBytes   = fread(fid,1,'uint64');
header.fileType        = char(fread(fid,3,'char')');
header.fileVersion     = fread(fid,1,'char');

if header.fileVersion < 3
    
    % event name of stream
    if header.fileVersion == 2
        header.eventName  = char(fread(fid,4,'char')');
    else
        header.eventName  = fliplr(char(fread(fid,4,'char')'));
    end
    
    % current channel of stream
    header.channelNum        = fread(fid, 1, 'uint16');
    % total number of channels in the stream
    header.totalNumChannels  = fread(fid, 1, 'uint16');
    % number of bytes per sample
    header.sampleWidthBytes  = fread(fid, 1, 'uint16');
    reserved                 = fread(fid, 1, 'uint16');
    
    % data format of stream in lower four bits
    header.dForm      = ALLOWED_FORMATS{bitand(fread(fid, 1, 'uint8'),7)+1};
    
    % used to compute actual sampling rate
    header.decimate   = fread(fid, 1, 'uint8');
    header.rate       = fread(fid, 1, 'uint16');
    
    % reserved tags
    reserved = fread(fid, 1, 'uint64');
    reserved = fread(fid, 2, 'uint16');
end

if header.fileVersion > 0
    % determine data sampling rate
    header.Fs = 2^(header.rate)*25000000/2^12/header.decimate;
end

fclose(fid);