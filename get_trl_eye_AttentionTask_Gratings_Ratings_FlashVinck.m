function [trl, taccept, trldesc] = get_trl_eye_AttentionTask_Gratings_Ratings_FlashVinck(cfg, trialEvents, behav)
%          startEvent    triggerEvent    endEvent
% --------------|-------------|-------------|---------------
%  preStartTime               0                  postEndTime

if ~exist('cfg', 'var'); cfg = []; end

% sampling rate of to be snippeted data
if ~isfield(cfg, 'samplingRate'); cfg.samplingRate = 1000; end

% trial defining events
if ~isfield(cfg, 'startEvent'); cfg.startEvent = 'TRIALID'; end
if ~isfield(cfg, 'endEvent'); cfg.endEvent = 'TRIAL_RESULT'; end
if ~isfield(cfg, 'triggerEvent'); cfg.triggerEvent = 'CONDITION'; end

% extra time before/after start/end event
if ~isfield(cfg, 'preStartTime'); cfg.preStartTime = 0; end
if ~isfield(cfg, 'postEndTime'); cfg.postEndTime = 0; end


if ~isfield(cfg, 'includedOutcomes')
    cfg.includedOutcomes = {'OK'};
end

% conversion factor from trialEvents to samples
convFac = 1000/cfg.samplingRate;

% get trial events with status 0 or 2
trialind =find(strncmp(trialEvents{3}, 'TRIAL_RES', 9) & ...
    (strncmp(trialEvents{4}, '2', 5) | ...
     strncmp(trialEvents{4}, '3', 5) | ...
     strncmp(trialEvents{4}, '4', 5) | ...
     strncmp(trialEvents{4}, '0', 5)));
endEvents = trialEvents{4}(trialind(2:end));
endEvents = cellfun(@(x) str2num(x), endEvents);

% event definitions
trialOutcome = num2cell(behav.trialOutcome);

% check if the eye data matches tdt events
assert(all(endEvents==[trialOutcome{:}]'),'eye trials dont match')
endInd = trialind(2:end);
nTrialsEdf = length(trialOutcome);

% putative start ind
startind =find(strncmp(trialEvents{3}, 'TRIALID', 9));
startInd = [];
for ii=1:nTrialsEdf
    diffInd = (startind-endInd(ii));
    [~, sInd] = max(diffInd(diffInd<0));
    startInd = [startInd; startind(sInd)];
end

n = 0;
for iTrialEdf = 1:nTrialsEdf
    n=n+1;
    %     if trialOutcome{ii} == 0
    % event MSG for current trial
    msg = trialEvents{3}(startInd(iTrialEdf):endInd(iTrialEdf));
    msg2 = trialEvents{4}(startInd(iTrialEdf):endInd(iTrialEdf));
    smpl = round(trialEvents{2}(startInd(iTrialEdf):endInd(iTrialEdf))); %sampling rate is 500 ===== /2
    
    if ~any(strncmp(msg, 'CONDITION', 10))
        if trialEvents{3}(endInd(iTrialEdf)) == 0; error('check this'); end
    end
    
    % find timings
    e.startEvent = smpl(strcmp(msg, cfg.startEvent));
    e.triggerEvent = smpl(strcmp(msg, cfg.triggerEvent));
    e.endEvent = smpl(strcmp(msg, cfg.endEvent));
    
    % start of trial
    trial(n).start = double(e.startEvent - cfg.preStartTime*cfg.samplingRate);
    
    % end of trial
    trial(n).end = double(e.endEvent + cfg.postEndTime*cfg.samplingRate);
    
    % trigger event within trial
    if isempty(e.triggerEvent)
        error('Trigger empty on trial %i', iBhvTrial)
    end
    
    if length(e.triggerEvent)>1
        trial(n).trigger = trial(n).start - e.triggerEvent(2);
    else
        trial(n).trigger = trial(n).start - e.triggerEvent(1);
    end
    
    % time of task change onset from cueon
    %     trial(n).stimOn = smpl(strcmp(msg, 'TRIALID')) + trial(n).trigger-trial(n).start;
    %     trial(n).respOn = smpl(strcmp(msg, 'FIXATION')) + trial(n).trigger-trial(n).start;
    trial(n).trialOutcome = trialOutcome{iTrialEdf};
    %     end
    
    % create trl matrixs
    trl = zeros(length(trial), 3);
    trl(:,1) = [trial.start]';
    trl(:,2) = [trial.end]';
    % trl(:,3) = [trial.trialOutcome]';
    trl(:,3) = [trial.trigger]';
    % trl(:,5) = [trial.stimOn]';
    % trl(:,6) = [trial.respOn]';
    
    % description of trl matrix columns
    trldesc = {
        'start'
        'end'
        %     'status'
        'trigger'
        %     'stimOn'
        %     'respOn'
        };
end
taccept = vertcat(trialOutcome{:});
end

function tokens = regexp_token_stripped(C,expr)
tokens = regexp(C,expr,'tokens');
tokens = cat(1,tokens{:});
tokens = cat(1,tokens{:});
end

