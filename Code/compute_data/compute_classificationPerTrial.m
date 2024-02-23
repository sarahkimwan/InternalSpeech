
%% compute classification per trial

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); %add folder to search path 
clc
clear all
%close all

%%

subjectId= 's2';
unit_region = 'SMG'; %all, PMV, SMG, AIP, or BA5, or S1X

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd,'Database', subjectId,'Data');

tableName = ['Table_' subjectId  '_Speech.mat'];

Data = load(fullfile(dataFolderPath,tableName)); %getting slighlty better results without ratefilt -> keep it?
flagShuffleTest = false; % if true, shuffle lables
flagSaveData = false; % saves the data 
flagLeaveOneOut = true; %if true = yes, if not true, then 8 fold cross validation
taskCue = 'All'; %All, Auditory, Written
PCA_variance_explained = 95; 

if flagShuffleTest
    numRep = 100; 
else
    numRep = 1; 
end
   
Go_data = Data.Go_data;    

if strcmp(taskCue, 'All')
    %keep all
elseif strcmp(taskCue, 'Auditory')
    Go_data = Go_data(ismember(Go_data.cueType, 'sound'),:);
elseif strcmp(taskCue, 'Written')
    Go_data = Go_data(ismember(Go_data.cueType, 'writing'),:);
end 

sessions = unique(Go_data.session_date);
numSessions = length(sessions);

flagPCA = true;
labels = Go_data.GoLabels;

numPhases = 6; 

numClasses = numel(unique(labels));
maxClassNumPerSessionTmp = zeros(numSessions,1);

%correctly predeclare variable size 
for n_session = 1:numSessions
    labels_tmp = labels(ismember(Go_data.session_date,sessions{n_session}));
    [counts, groupnames] = groupcounts(labels_tmp);
    maxClassNumPerSessionTmp(n_session) = max(counts);
end

maxClassNumPerSession = max(maxClassNumPerSessionTmp);

if flagLeaveOneOut
    k = maxClassNumPerSession*numClasses;  
else
    k = 8;
end

%predeclare variables 

errTest = nan*ones(numRep,k,numPhases,numSessions); 
keepTestingLabels = nan*ones(numRep,maxClassNumPerSession*numClasses,numPhases,numSessions);  
keepPredictedLabels = nan*ones(numRep,maxClassNumPerSession*numClasses,numPhases,numSessions); 

PhaseNames = {'ITI', 'Cue','Delay1','ImaginedSpeech', 'Delay2', 'Speech'};

if strcmp(unit_region, 'S1X') && strcmp(subjectId, 's3')
    sessionToInclude = setdiff(1:numSessions,2); %no data for session 2 in S1
else
    sessionToInclude = 1:numSessions; 
end


for n_session =  sessionToInclude
    disp(['Classifying session ' sessions{n_session} ]);
    idxThisSession = ismember(Go_data.session_date, sessions(n_session));
    
    if strcmp('SMG', unit_region)
        SessionData = Go_data.SMG_Go(idxThisSession,:);
    elseif strcmp('S1X', unit_region)
        SessionData = Go_data.S1X_Go(idxThisSession,:);
    end  
    
    labelsPerSession = labels(idxThisSession);  
    timePhaseLabels = Go_data.time_phase_labels(idxThisSession); 
       
    if numSessions ==1
        % Take all session days together
        if strcmp('SMG', unit_region)
            SessionData = Go_data.SMG_Go(ismember(Go_data.session_date,sessions),:);
        
        elseif strcmp('S1X', unit_region)
            SessionData = Go_data.S1X_Go(ismember(Go_data.session_date,sessions),:);
        end 
        labelsPerSession = labels(ismember(Go_data.session_date,sessions),:);
        timePhaseLabels = Go_data.time_phase_labels(ismember(Go_data.session_date,sessions),:);
    end
    

    for n_phase = 1:numPhases

        DataPerPhase = cell2mat(arrayfun(@(x,y) nanmean(x{1,1}(y{:}== n_phase,:),1),SessionData,timePhaseLabels, 'UniformOutput', false));
               
       [errTrainPhase,errTestPhase, labelsTest,predictedTest]= classification.LDA_classification_rep(DataPerPhase,labelsPerSession,...
         'flagPCA', true, 'PCA_variance', PCA_variance_explained, 'classifierType', 'diaglinear', 'numRep', numRep, 'flagLeaveOneOut',flagLeaveOneOut, ...
         'k', k, 'flagRandomPerm', flagShuffleTest, 'flagTuning', false,'flagErrorMatrix', false); 
        
       errTest(:,1:size(errTestPhase,2),n_phase,n_session) = errTestPhase;
       
       keepTestingLabels(:,1:size(labelsTest,2),n_phase,n_session) = labelsTest;
       keepPredictedLabels(:,1:size(predictedTest,2),n_phase,n_session) = predictedTest;
   
    end          
end 

keepMeanAcc = (1-squeeze(nanmean(errTest(1,:,:,:),2)))*100;

if flagSaveData
    
   dataFolderPath = fullfile(pwd,'Database',subjectId,'Data','ClassificationData','ClassificationPerTrial');

    if flagShuffleTest
        save(fullfile(dataFolderPath,[unit_region '_' taskCue  '_Speech_Shuffle']),'errTest', 'keepTestingLabels', 'keepPredictedLabels');

    else        
        save(fullfile(dataFolderPath,[unit_region '_' taskCue  '_Speech_Data']),'errTest', 'keepTestingLabels', 'keepPredictedLabels');
    end 
end 

%% quick plot

if ndims(keepMeanAcc) >2
    keepMeanAcc = squeeze(mean(keepMeanAcc,1));
end 
figure()

err_ci =utile.calculateErrCi(keepMeanAcc);

plot( squeeze(keepMeanAcc), 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none', 'Color', [0,0,0]);
hold on
hold on 

ColorModel = [62,164,144]/255;
 errorbar(1:size(keepMeanAcc,1), nanmean(keepMeanAcc,2), err_ci(1,:), err_ci(2,:),'-s', 'LineStyle', '--','MarkerSize',10,...
            'MarkerEdgeColor',ColorModel,'MarkerFaceColor',ColorModel , 'Color', ColorModel) 
      
xticks(1:numPhases);
xticklabels(PhaseNames)
xtickangle(-45)

ylabel('Classification accuracy [%]');
title(['Classification using ' taskCue ' for training - ' unit_region]);
ylim([0 100])   
meanAll = mean(keepMeanAcc,2);
for n_phsdr = 1:numPhases
    text(n_phsdr+0.2, meanAll(n_phsdr),[num2str(meanAll(n_phsdr),3) '%'])
end 
