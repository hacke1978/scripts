function dataSnip = convert_edf_fieldtrip(cfg,file)
addpath('/opt/fieldtrip_esi');
% this function makes 3 txt files that can later be read into Matlab
% using the matlab function dlmread(), the files are space ' ' delimited
% when there are blinks, eyelink puts in '.'
% this function replaces them with Nan
% it is up to the user to decide what to do with the Nans

% grep segments originally by Jarrod Dowell
% rewritten by Katharine Shapcott to convert straight to fieldtrip format 06/12/2016

if ~exist('file', 'var')
    [filename, pathname] = uigetfile(...
        {'*.asc'}, 'Pick a converted edf file', '~/confidence_psychophysics/edf');
    if ~ischar(filename); return; end % exit if cancel
    file = fullfile(pathname, filename);
    clear filename pathname
end

if ~exist('cfg','var'); cfg = []; end
if ~isfield(cfg,'savepath');
    [cfg.savepath, ~] = fileparts(file);
end

% get file parts
[pathstr,name,ext] = fileparts(file);
subStr = vertcat(strsplit(name, '_'));
[subNum, ~, ~, subDate, ~] = deal(subStr{:});
% subDate = name(end-12:end);
% subNum = regexp(name,'Sub(\d+)\d\d-', 'tokens'); %only way to parse out participant number, using the start of the date portion
% subNum = str2double(subNum{1}{1});

fprintf('Converting edf file %s to fieldtrip format\n', name)

% add '' to file incase there are any spaces, grep will thank you
file = ['''',file,''''];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       read and export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MSG events
fprintf('Read MSG events\n')

if ismac
    [s,trlid] = system(['grep "MSG" ' [file] ' | awk ''{ print $1,$2,$3; }'' | sed -e s/$//g']);
elseif ispc
    [s,trlid] = system(['grep! "MSG" ' [file] ' | awk ''{ print $1,$2,$3; }'' | sed -e s/$//g']);
else
    [s,trlid] = system(['grep "MSG" ' [file] ' | awk ''{ print $1,$2,$3,$4,$5; }'' | sed -e s/$//g']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         write file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not slow at all
%fid = fopen([pathstr,'/',name,'.trledf'], 'w'); % trial number should match corresponding bhv
%fprintf(fid,'%s',trlid);
%fclose('all');

tmpTrialId = textscan(trlid,'%s%d%s%s%s');


%% Continuous eye data

fprintf('Read continuous eye data\n')

% read anything that starts with two numbers seperated by whitespace
if ismac
    [s,cont] = system(['grep ''^.[0-9][0-9]*[[:space:]][0-9]'' ' [file] ' | awk ''{ print $1,$2,$3,$4,$5; }'' | sed -e s/$//g']);
else
    [s,cont] = system(['grep ''^.[0-9][0-9]*[[:space:]][0-9]'' ' [file] ' | awk ''{ print $1,$2,$3,$4,$5; }'' | sed -e s/$//g']);
end

% find any '.' and replace with 'Nan'
% is any ' . ', these indicies will all be shifted because of the
% spaces, the spaces distinguish '9.8' from ' . ' (i.e., decimal places)
% thus add 1 to account for the shift

% because I am adding 'Nan', then each char in the whole string
% will be shifted by and additional 2 positions relative to its neighbour
eyeblanks = strfind(cont,' . ');


if ~isempty(eyeblanks)
    nBlanks = length(eyeblanks);
    % let the user know, they can deal with the Nans as they wish
    disp(['Found missing data...(',num2str(nBlanks),' values)... replacing with Nans']);
    
    cont(eyeblanks+1)='';   %remove '.' from the data
    
    eyeblanks = eyeblanks+(0:2:(nBlanks*2)-2)+1;
    
    % now that I have accounted for the shifts I can loop a string insert
    % algorthim
    
    % KS- loop was too slow, use 'blanks' to create new vector and then
    % fill in with 'N''a''n' and the old data (without '.')
    
    nFinal = numel(cont)+numel(eyeblanks)*3;
    newstr = blanks(nFinal);
    
    newstr(eyeblanks) = 'N';
    newstr(eyeblanks+1) = 'a';
    newstr(eyeblanks+2) = 'n';
    
    nanindices = [eyeblanks eyeblanks+1 eyeblanks+2];
    newstr(setdiff(1:nFinal,nanindices)) = cont(1,:);
    cont = newstr;
    clear newstr
    
    %         for eaBl = 1:nBlanks
    %             ind = blanks(eaBl);
    %             tmp = [cont(1,1:ind-1),'Nan',cont(1,ind+1:end)];
    %             cont = tmp;
    %         end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         write file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not slow at all
%fid = fopen([pathstr,'/',name,'.contedf'], 'w'); % continous data
%fprintf(fid,'%s',cont);
%fclose('all');

%%

fprintf('Create eye data structure\n')

tmpCont = textscan(cont,'%d%f%f%f%s');


%% read header
hdr = [];
EDF = read_edf_header(file(2:end-1));

hdr.orig = EDF;

hdr.subjectNumber = subNum;
hdr.subjectDate = subDate;
hdr.trialEvents = tmpTrialId;
hdr.TimeStampPerSample = diff(tmpCont{1}([1,2]));
hdr.Fs = double(1000 / hdr.TimeStampPerSample);
hdr.nChans = 3;
hdr.FirstTimeStamp = double(tmpCont{1}(1));
hdr.nSamples = (tmpCont{1}(end) - tmpCont{1}(1) + 1)/hdr.TimeStampPerSample;
hdr.nSamplesPre = 0;
hdr.nTrials = 1;


% fid=fopen(file(2:end-1));
% for iLine = 1:24
%     hdr{iLine,1} = fgetl(fid);
% end
% fclose(fid);

%% fill in the blanks

index = ((tmpCont{1} - hdr.FirstTimeStamp)/hdr.TimeStampPerSample) + 1;

dat = nan(hdr.nChans, hdr.nSamples);
dat(:,index) = cat(2,tmpCont{2:4})';

%% convert events to ms

hdr.trialEvents{2} = (hdr.trialEvents{2} - hdr.FirstTimeStamp) + 1;


%% calibrate data
% code from fieldtrip- we don't have PhysMax etc for some reason


% hdr.orig.Cal = (hdr.orig.PhysMax-hdr.orig.PhysMin)./(hdr.orig.DigMax-hdr.orig.DigMin);
% % Calibrate the data

% calib = diag(hdr.orig.Cal);
% dat = sparse(calib) * dat;

%% make raw data

data = [];
data.hdr = hdr;
data.label = {'eye_x', 'eye_y', 'pupil'};
data.time = {(0:(hdr.nSamples-1))/hdr.Fs}; %(tmpCont{1}(1):tmpCont{1}(end))-tmpCont{1}(1)
data.trial = {dat};

data.fsample = hdr.Fs;
data.sampleinfo = [1 length(data.time{1})];

%% save data
fprintf('Save eye\n')

saveName = [cfg.savepath, '.eye'];
try
    save(saveName, 'data', '-v6')
catch me
    warning('Using -v7.3 flag for saving due to data length')
    save(saveName, 'data', '-v7.3')
end

%% snippet data

fprintf('Snippet data\n')

% load behavioural

% behav = load([file(2:end-5) '.mat']);
%conditions = behav.cfg.lookup(behav.cfg.recording(1:behav.cfg.validTrials,2),:);


% calculate trial starts

trlCfg = [];
trlCfg.samplingRate = data.fsample;
redefCfg = [];
[redefCfg.trl, trldesc] = get_trl_eye_AttentionTask_Gratings_Ratings_FlashVinck(trlCfg, data.hdr.trialEvents, cfg);

dataSnip = ft_redefinetrial(redefCfg, data);
dataSnip.taccept = trldesc;
dataSnip.hdr = hdr;
dataSnip.hdr.label = data.label;

% %% add screen info
% 
% dataSnip.cfg.behavInfo.fixWindow = behav.cfg.fixWindow;
% dataSnip.cfg.behavInfo.screenRect = behav.cfg.screenRect;
% dataSnip.cfg.behavInfo.screenHeight = behav.cfg.screenHeight;
% dataSnip.cfg.behavInfo.screenDistance = behav.cfg.screenDistance;

% %% save snippetted data
% fprintf('Save stimOn.eye\n')
% 
% saveName = fullfile(cfg.savepath, [name, '.stimOn.eye']);
% try
%     save(saveName, 'data', '-v6')
% catch me
%     warning('Using -v7.3 flag for saving due to data length')
%     save(saveName, 'data', '-v7.3')
% end

% %% test plots
% testPlot = 0;
% 
% if testPlot
%     h = figure;
%     for k = 1:length(dataSnip.trial)
% %         bhvTrl = dataSnip.trialinfo(k,1);
%         if dataSnip.taccept(k) == 0
%             plot(dataSnip.trial{k}(1,:),dataSnip.trial{k}(2,:), 'ko')
%         else 
%             plot(dataSnip.trial{k}(1,:),dataSnip.trial{k}(2,:), 'r.')
%             end
%         hold on
%         %plot(bhv.AnalogData{bhvTrl}.EyeSignal(round(bhv.CodeTimes{bhvTrl}(1)/2):round(bhv.CodeTimes{bhvTrl}(end)/2),:) *100)
% %         plot([sum(dataSnip.trialinfo(k,9)),sum(dataSnip.trialinfo(k,9))], ylim, 'k')
% %         waitfor(h)
%     end
% end



end


%% subfunctions

function EDF = read_edf_header(filename)

fid=fopen(filename);
for iLine = 1:25
    fLine{iLine,1} = fgetl(fid);
end
fclose(fid);

EDF.FILE.FID = fid;
EDF.FILE.OPEN = 1;
[EDF.FILE.Path, EDF.FILE.Name, EDF.FILE.Ext] = fileparts(filename);
EDF.FileName = filename;

EDF.CONVERTED=fLine{1};
EDF.DATE=fLine{2};
EDF.TYPE=fLine{3};
EDF.VERSION=fLine{4};
EDF.SOURCE = fLine{5};
EDF.FULLVERSION=fLine{6};
EDF.CAMERA=fLine{7};
EDF.SERIALNUMBER=fLine{8};
EDF.CAMERA_CONFIG=fLine{9};

EDF.RECCFG = fLine{13};
EDF.ELCLCFG = fLine{14};
EDF.GAZE_COORDS = fLine{15};
EDF.THRESHOLDS = fLine{16};
EDF.ELCL_PROC = fLine{17};
EDF.ELCL_PCR_PARAM = fLine{18};
EDF.MODE = fLine{19};
EDF.START = fLine{20};
EDF.PRESCALER = fLine{21};
EDF.VPRESCALER = fLine{22};
EDF.PUPIL = fLine{23};
EDF.EVENTS = fLine{24};
EDF.SAMPLES = fLine{25};


end

