
%% Compute tuned units per phase for each session day

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

% Add current folder and its subfolders to the search path
addpath(genpath(pwd)); 

clc
clear all 
%close all

%%

subjectId = 's2'; %s3
unitRegion = 'SMG'; %S1X

flagRegressionTuning = true; %true = compute linear regression analysis. False = compute kruskalwallis test

% Define the data folder path based on subject ID
saveDataLocation = fullfile(pwd,'Database', subjectId,'Data');

tableName = ['Table_' subjectId  '_Speech.mat'];
 
flagBinPerBin = false;
multipleComparePhase = true; 

flagSaveData = false; 

%cue type per participant
if strcmp(subjectId, 's2')
   taskCuesAll = {'Auditory', 'Written'};
elseif strcmp(subjectId, 's3')
    taskCuesAll = {'Written'}; %only Written cue trials for s3
else
    keyboard % unknown subject
end

Data = load(fullfile(saveDataLocation,tableName));

Sessions = unique(Data.Go_data.session_date);
numSessions = length(Sessions);

%predeclare variables 
numCues = numel(taskCuesAll);
tuned_channels_per_phase = cell(numCues, numSessions);
tuned_channels_per_phase_vector = cell(numCues, numSessions);
sum_bin_all = cell(numCues, numSessions);
numUnitsPerSession = zeros(numCues,numSessions);
timePhaseLabelsPerSession = zeros(numCues,numSessions,157);

for n_cue = 1:numCues
    
    taskCue = taskCuesAll{n_cue};
    Go_data = Data.Go_data;    

    if strcmp(taskCue, 'All')
        %keep all
    elseif strcmp(taskCue, 'Auditory')
        Go_data = Go_data(ismember(Go_data.cueType, 'sound'),:);
    elseif strcmp(taskCue, 'Written')
        Go_data = Go_data(ismember(Go_data.cueType, 'writing'),:);
    end

    labels = Go_data.GoLabels;

    for n_session = 1:numSessions

        disp(['Classifying session ' Sessions{n_session} ]);
        idxThisSession = ismember(Go_data.session_date, Sessions(n_session));

        if strcmp('SMG', unitRegion)
            sessionData = Go_data.SMG_Go(idxThisSession,:);
        elseif strcmp('S1X', unitRegion)
            sessionData = Go_data.S1X_Go(idxThisSession,:);
        end
        
        labelsPerSession = labels(idxThisSession);
        timePhaseLabels = Go_data.time_phase_labels(idxThisSession);
        numUnitsPerSession(n_cue,n_session) = size(sessionData{1},2);
        timePhaseLabelsPerSession(n_cue, n_session,:) = timePhaseLabels{1};

       %Compute index of units that are tuned
       if flagRegressionTuning
           
             [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin,numTunedChannelsPerCategory,tunedChannelsPerPhasePerCategoryTmp,...
                 ~,~,p_per_phase] ...
                  = classification.getRegressionTunedChannels(sessionData,labelsPerSession, ...
                     timePhaseLabels, 'multcompare', multipleComparePhase, 'BinperBinTuning', flagBinPerBin);

                 tunedChannelsPerPhasePerCategory{n_cue,n_session} = tunedChannelsPerPhasePerCategoryTmp;
                 
           condToTest = arrayfun(@(x) utile.image2class_simple(x),  unique(labelsPerSession), 'UniformOutput', false);

            tuned_channels_per_graps{n_cue,n_session} = numTunedChannelsPerCategory;
       else

           tuned_channels_per_graps{n_cue,n_session} = []; % not computed for Kurkalwallis test
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
    
    saveFolder = fullfile(pwd,subjectId,'TuningData');
    save(fullfile(saveFolder,[unitRegion '_' reg_sub{flagRegressionTuning+1}  '_Speech_tuning']),'sum_bin_all','tuned_channels_per_phase','tuned_channels_per_phase_vector','timePhaseLabels','numUnitsPerSession','tuned_channels_per_graps', ...
        'timePhaseLabelsPerSession', 'numSessions', 'numCues','labelsPerSession','tunedChannelsPerPhasePerCategory');
end 
