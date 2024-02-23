%% Calculate cross decoding analysis Audio - Written cue 

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

clc
clear all 
%close all

%% 
unit_region = 'SMG'; %SMG or S1X
subjectId = 's2';
trainingCue = 'Written'; %switch between Auditory and Written
%trainingCue = 'Auditory';

% Define the data folder path based on subject ID
saveDataLocation = fullfile(pwd,'Database', subjectId,'Data');

tableName = ['Table_' subjectId  '_Speech.mat'];

flagSaveData = false; 
flagShuffleTest = false ; 

if flagShuffleTest
    numReps = 100; 
    %repeat shuffle distribution 100 times to generate null distribution
else
    numReps = 1; 
end

if strcmp(trainingCue, 'Written')
    trainingCueModality = 'writing'; %written ,speech
    testingCueModality = 'sound';
    testingCue = 'Auditory';
    ColorsModel = [30,136,229]/255;
else
    trainingCueModality = 'sound'; %written ,speech
    testingCueModality = 'writing';
    testingCue = 'Written';
    ColorsModel = [93,191,59]/255;

end 

%load data
Data = load(fullfile(saveDataLocation,tableName));

trainingData = Data.Go_data(ismember(Data.Go_data.cueType, trainingCueModality),:);
testingData = Data.Go_data(ismember(Data.Go_data.cueType, testingCueModality),:);
Sessions = unique(trainingData.session_date);
numSessions = length(Sessions);

flagPCA = true;
trainingLabels = trainingData.GoLabels;
testingLabels = testingData.GoLabels;
classifierType = 'diaglinear'; 

numPhases = 6; 

%predeclare variables 
errTrain = nan*ones(numReps, numPhases,numSessions); 
errTest = nan*ones(numReps, numPhases,numSessions); 

keepTestingLabels = cell(numReps,numPhases, numSessions);
keepPredictedLabels = cell(numReps, numPhases,numSessions);

%save number of PC per session
keepNumPCAUnits = nan*ones(numSessions,1);
 
phaseNames = {'ITI', 'Cue','Delay1','ImaginedSpeech', 'Delay2', 'Speech'};


%compute cross-classification for each session day
for n_session = 1:numSessions
    disp(['Classifying session ' Sessions{n_session} ]);
    
    idxTrainingSession = ismember(trainingData.session_date, Sessions(n_session));
    idxTestingSession = ismember(testingData.session_date, Sessions(n_session));

    %extract data of session day
    if strcmp('SMG', unit_region)
        sessionDataTrain = trainingData.SMG_Go(idxTrainingSession,:);
        sessionDataTest = testingData.SMG_Go(idxTestingSession,:);

    elseif strcmp('S1X', unit_region)
        sessionDataTrain = trainingData.S1X_Go(idxTrainingSession,:);
        sessionDataTest = testingData.S1X_Go(idxTestingSession,:);
    end    
        
    labelsTrain = trainingLabels(idxTrainingSession);
    labelsTest = testingLabels(idxTestingSession);
     
    timePhaseLabelsTrain = trainingData.time_phase_labels(idxTrainingSession);
    timePhaseLabelsTest = testingData.time_phase_labels(idxTestingSession); 

    for n_rep  = 1:numReps
        
        %randomize training labels to compute shuffle distribution
        if flagShuffleTest
            labelsTrain = labelsTrain(randperm(length(labelsTrain)));
        end 

        %compute classification for each of six task phases
        for n_phase = 1:numPhases
            
            dataTrain = cell2mat(arrayfun(@(x,y) nanmean(x{1,1}(y{:}== n_phase,:),1),sessionDataTrain,timePhaseLabelsTrain, 'UniformOutput', false));
            dataTest = cell2mat(arrayfun(@(x,y) nanmean(x{1,1}(y{:}== n_phase,:),1),sessionDataTest,timePhaseLabelsTest, 'UniformOutput', false));
           
            [errTrain(n_rep,n_phase,n_session),errTest(n_rep,n_phase,n_session), ~, predictedTest] = ...
                classification.LDA_classification_train_test(dataTrain, dataTest,labelsTrain, labelsTest, ...
              'flagPCA', flagPCA, 'PCA_variance', 90, 'classifierType', classifierType, 'flagTuning', false); 

            keepTestingLabels{n_rep,n_phase, n_session} = labelsTest;
            keepPredictedLabels{n_rep,n_phase, n_session} = predictedTest;                    
  
        end 
     end                
end

keepLabels = struct; 
keepLabels.keepTestingLabels = keepTestingLabels; 
keepLabels.keepPredictedLabels = keepPredictedLabels;

plotData = struct;
plotData.TrainingCue = trainingCue; 
plotData.TestingCue = testingCue; 
plotData.numPhases = numPhases; 
plotData.ColorsModel = ColorsModel;
plotData.phaseNames = {'ITI', 'Cue', 'D1', 'Imagined', 'D2', 'Speech'};
plotData.Sessions = Sessions; 

if flagSaveData

      saveDataLocation = fullfile(pwd,'Database','s2','Data','ClassificationData','AudioWrittenCrossDecoding');
    
    if flagShuffleTest
        save(fullfile( saveDataLocation, [unit_region '_' trainingCue 'TrainingData_Shuffle' ]), 'errTest','plotData', 'keepLabels');
    else
        save(fullfile( saveDataLocation, [unit_region '_' trainingCue 'TrainingData' ]), 'errTest','plotData', 'keepLabels');
    end 
end

%% plot results (quick plot)

if ndims(squeeze(errTest)) > 2
    errTest = squeeze(mean(errTest,1));
end 
    
x = 1:numPhases;                                  
acc = (1 - squeeze(errTest)')*100;
NA = size(acc,1);
meanAcc = mean(acc);
keepMeanAcc = meanAcc;
SEM_A = std(acc) / sqrt(NA); % Standard Error Of The Mean
CI95_A = SEM_A * tinv(0.975, NA-1);  
keepCI = CI95_A;

pretend_online = figure('units','normalized','outerposition',[0 0 0.2 0.4]);
plot( acc', 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none', 'Color', [0 0 0]);
hold on
e = errorbar(1:length(meanAcc), meanAcc, CI95_A,'-s' ,'LineStyle', '--','MarkerSize',10,...
            'MarkerEdgeColor',ColorsModel,'MarkerFaceColor',ColorsModel , 'Color', ColorsModel) ;
ylabel('Classification accuracy [%]')
title(['Training data : ' trainingCue]);
legend(e ,['Testing data : ' testingCue]);

xticks(1:numPhases)
xticklabels(phaseNames)
xtickangle(-45)
xlim([0.5, 6.5])
set(gca, 'fontSize', 12)

