
%% Compute tuned units per 50ms timebins for each session day - aligned to audio onset

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway
% Add current folder and its subfolders to the search path
addpath(genpath(pwd)); 

clc
clear all 
%close all
%%

subjectId = 's3';
unitRegion = 'SMG'; 

dataFolderPath = fullfile(pwd, 'Database', subjectId, 'Data');

flagRegressionTuning = false; 

tableName = ['Table_' subjectId '_Speech.mat'];
 
flagBinPerBin = true;
multipleComparePhase = true; 
flagTunedChannels = true;

flagSaveData = false; 

%chose cue type:
if strcmp(subjectId, 's2')
   taskCuesAll = {'Auditory', 'Written'};
elseif strcmp(subjectId, 's3')
    taskCuesAll = {'Written'}; %only Written cue trials for s3
else
    keyboard % unknown subject
end

numCues = numel(taskCuesAll);
Data = load(fullfile(dataFolderPath,tableName));
Sessions = unique(Data.Go_data.session_date);
numSessions = length(Sessions);

%predeclare variables 
tuned_channels_per_phase = cell(numCues, numSessions);
tuned_channels_per_phase_vector = cell(numCues, numSessions);
sum_bin_all = cell(numCues, numSessions);
numUnitsPerSession = zeros(numCues,numSessions);
timePhaseLabelsPerSession = zeros(numCues,numSessions,157);

audio_start = [47	48	47	51	49	49	48	48	48	49	49	50	49	48	50	48	49	48	50	52	...
    48	48	50	48	48	48	49	49	48	48	51	48	48	48	50	49	48	49	51	48	49	47	50	48	...
    51	50	48	49	49	49	50	47	50	47	50	48	48	48	47	50	48	47	50	48	48	49	48	47	50 ...
    48	48	48	48	47	47	47	47	49	49	48	47	51	51	49	47	48	47	47	49	48	47	48	47	47	...
    47	49	48	47	51	48	51	46	49	47	49	48	48	47	48	49	47	48	47	47	48	48	48	48	48	...
    49	47	48	49	48	50	50	48	49	48	48	47	47	48	48	49	49	49	49	49	48	48	51	48	48	51 ...
    47	47	50	48	51	49	48	48	49	50	48	48	48	49	50	48	47	49	48	47	48	47	47	47	48	48 ...
    48	50	49	49	49	48	48	49	50	50	48	47	48	48	48	50	50	47	55	54	54	54	54	54	54	55 ...
    55	54	55	54	55	54	54	55	54	54	54	54	54	54	54	54	55	54	53	54	54	54	55	55	54	55 ...
    55	54	55	54	55	54	55	54	54	55	54	55	53	54	55	55	55	55	54	55	54	54	54	54	54	54	...
    54	54	54	55	55	54	54	53	54	54	54	53	55	54	53	54	54	54	54	54	54	53	54	54	53	54	...
    54	54	54	54	54	55	53	54	55	54	55	53	55	54	54	53	54	53	53	55	55	54	54	54	54	55	...
    54	55	54	54	53	55	54	54	53	53	54	54	54	54	54	54	52	54	54	55	54	55	54	54	54	55	...
    54	54	53	55	53	54	55	54	54	54	54	54	54	54	54	54	54	55	54	55	54	55	54	54	54	54	...
    55	54	54	55	53	54	54	53	54	54	54	54	54	54	54	55	51	55	54	54	54	53	54	54	54	54	...
    54	54	53	54	54	54	54	54	54	54	54	54	54	54	55	54	53	55	55	54	54	54	53	53	54	54	...
    54	55	54	54	53	54	55	54	55	53	54	54	55	54	54	55	55	54	55	54	54	53	53	54	53	54	...
    53	53	54	54	54	55	54	54	54	55	54	54	54	54	53	54	54	53	54	54	55	55	54	53	54	54	...
    54	55	54	54	54	53	54	55	54	53	54	55	54	54	54	54	54	54	55	53	54	55	54	53	53	53 54	...
    54	54	54	55	54	54	55	55	54	55	55	54	55	55	55	55	55	54	53	54	55	54	53	54	55  54	52 ...
    53	54	54	54	55	54	55	54	55	55	53	54	54	53	53	55	54	53	54	53	55	54	54	55	54 54 54 ...
    54	54	54	54	54	55	54	54	54	54	54	54	54	54	54	55	54	54	54	55	54	54	54	54	54	55	54	54  ];

written_start(1:size(Data.Go_data,1)) = 42;
maxLength = 18; 
audio_start(end:size(Data.Go_data,1)) = 54;

%loop over cue conditions  

for n_cue = 1:numCues
    
    taskCue = taskCuesAll{n_cue};
    Go_data = Data.Go_data; 
    time_phase_labels_fixed = Go_data.time_phase_labels{1};
   
   if strcmp(taskCue, 'All')
        %keep all
    elseif strcmp(taskCue, 'Auditory')
        Go_data = Go_data(ismember(Go_data.cueType, 'sound'),:);
    elseif strcmp(taskCue, 'Written')
        Go_data = Go_data(ismember(Go_data.cueType, 'writing'),:);
   end
  
    labels = Go_data.GoLabels;

    for n_session =  1:numSessions
        disp(['Classifying session ' Sessions{n_session} ]);
        idxThisSession = ismember(Go_data.session_date, Sessions(n_session));

        if strcmp('SMG', unitRegion)
            sessionData = Go_data.SMG_Go(idxThisSession,:);
        elseif strcmp('S1X', unitRegion)
            sessionData = Go_data.S1X_Go(idxThisSession,:);
        end
        
        labelsPerSession = labels(idxThisSession);
        timePhaseLabels = cell(nnz(labelsPerSession), 1);
        timePhaseLabels(:) = {time_phase_labels_fixed};
        flagSelectWords = false;

        
        if strcmp(taskCue, 'Auditory')
            %align audio neural data to audio cue start
            cueStartPerSession = audio_start(idxThisSession);
            cuePhaseEndIdx = cueStartPerSession + maxLength;
            cuePhaseIdxAll{n_cue,n_session}= arrayfun(@(x,y) x:y, cueStartPerSession, cuePhaseEndIdx, 'UniformOutput',false); 
        end

        numUnitsPerSession(n_cue,n_session) = size(sessionData{1},2);
        timePhaseLabelsPerSession(n_cue, n_session,:) = timePhaseLabels{1};
       
        if strcmp(taskCue, 'Auditory')
        
            % Align auditory neural data to cue onset
            cueIdxTmp = cellfun(@(x) find(x == 2), timePhaseLabels, 'UniformOutput',false);
            
            notToKeep = cellfun(@(x,y) ~ismember(x,y), cueIdxTmp,cuePhaseIdxAll{n_cue,n_session}', 'UniformOutput',false);
            idxToDelete = cellfun(@(x,y) x(y), cueIdxTmp,notToKeep, 'UniformOutput',false);
            idxToKeep = cellfun(@(x,y) setdiff(1:size(x,1), y), sessionData,idxToDelete, 'UniformOutput',false);
            diffLength = size(cueIdxTmp{1},1) - size(idxToDelete{1});
            sessionDataTmp =  cellfun(@(x,y) x(y,:), sessionData,idxToKeep, 'UniformOutput',false);
            % Replace shifted values by nan
            sessionData = cellfun(@(x) [x(1:cueIdxTmp{1}(1)+maxLength,:) ; nan*zeros(size(idxToDelete{1},1), size(x,2)); x(cueIdxTmp{1}(1)+maxLength+1:end,:)], sessionDataTmp, 'UniformOutput',false  );
        
        end 


           %Compute index of units that are tuned
           if flagRegressionTuning   

                [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin,numTunedChannelsPerCategory, ...
                    tunedChannelsPerPhasePerCategoryTmp,~,~,p_per_phase] ...
                = classification.getRegressionTunedChannels(sessionData,labelsPerSession, ...
                 timePhaseLabels, 'multcompare', multipleComparePhase, 'BinperBinTuning', flagBinPerBin);
                
                tunedChannelsPerPhasePerCategory{n_cue,n_session} = tunedChannelsPerPhasePerCategoryTmp;
                tuned_channels_per_graps{n_cue,n_session} = numTunedChannelsPerCategory;

           else
                tuned_channels_per_graps{n_cue,n_session} = [];
                [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin]= classification.getTunedChannels(sessionData,labelsPerSession, ...
                timePhaseLabels, 'multcompare', multipleComparePhase,'removeITItuning', 'false', 'BinperBinTuning', flagBinPerBin);   

                sumBin = sumBin';
           end 
           
           if nnz(sumBin) > 0
               sum_bin_all{n_cue, n_session } = sumBin;
           else
               sum_bin_all{n_cue, n_session } = [];
           end

            tuned_channels_per_phase{n_cue,n_session} = sumPhase;
            tuned_channels_per_phase_vector{n_cue,n_session} = tunedChannelsPhase;
    end  
end 

if flagSaveData
   
    reg_sub = {'KruskalWallis','Regression'};
    
    saveFolder = fullfile(pwd,subjectId,'AlignedTuningData');

    if flagRegressionTuning
         save(fullfile(saveFolder,[unitRegion '_' reg_sub{flagRegressionTuning+1}  '_Speech_tuning']),'sum_bin_all','tuned_channels_per_phase','tuned_channels_per_phase_vector','timePhaseLabels','numUnitsPerSession','tuned_channels_per_graps', ...
        'timePhaseLabelsPerSession', 'numSessions', 'numCues','labelsPerSession','tunedChannelsPerPhasePerCategory','audio_start', 'written_start','cuePhaseIdxAll');

    else
        save(fullfile(saveFolder,[unitRegion '_' reg_sub{flagRegressionTuning+1}  '_Speech_tuning']),'sum_bin_all','tuned_channels_per_phase','tuned_channels_per_phase_vector','timePhaseLabels','numUnitsPerSession','tuned_channels_per_graps', ...
        'timePhaseLabelsPerSession', 'numSessions', 'numCues','labelsPerSession','audio_start', 'written_start','cuePhaseIdxAll');

    end
   
end 
