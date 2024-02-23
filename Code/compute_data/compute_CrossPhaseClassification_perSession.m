

%% compute cross-classification per trial

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); %add folder to search path 
clc
clear all
%close all

%%

subjectId= 's2';
unit_region = 'SMG'; 

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd,'Database', subjectId,'Data');

tableName = ['Table_' subjectId  '_Speech.mat'];

Data = load(fullfile(dataFolderPath,tableName)); 

taskCue = 'Written'; %Auditory

Go_data = Data.Go_data;

if strcmp(taskCue, 'All')
    %keep all
elseif strcmp(taskCue, 'Auditory')
    Go_data = Go_data(ismember(Go_data.cueType, 'sound'),:);
elseif strcmp(taskCue, 'Written')
    Go_data = Go_data(ismember(Go_data.cueType, 'writing'),:);
end 

sessions = unique(Go_data.session_date);
numSessions = numel(sessions);
k = 8;
classifier = 'diaglinear'; 
numRepetition = 1; 
numPhases = 6; 

%predeclare variables 
flagLeaveOneOut = true; 
PCA_variance =  95;
flagPCA = true;
flagSaveData = false;

for n_session = 1:numSessions
    disp(['Classifying session ' sessions{n_session} ]);
    
    if numSessions ~= 1
        idxThisSession = ismember(Go_data.session_date, sessions(n_session));  
    else  
        idxThisSession = ismember(Go_data.session_date, sessions);

    end 
    
    sessionDataAll = Go_data(idxThisSession,:);
    labels = Go_data.GoLabels(idxThisSession);
    timePhaseLabelsAll = Go_data.time_phase_labels(idxThisSession);

    if strcmp('SMG', unit_region)
        SessionData= sessionDataAll.SMG_Go;
    elseif strcmp('S1X', unit_region)
        SessionData = sessionDataAll.S1X_Go;
    end
    
    timePhaseLabels =timePhaseLabelsAll{1};
    labelsPerSession = labels; 

    for rep = 1:numRepetition 

        dataPerPhaseAll = cell(numPhases,1);

        for n_testPhase = 1:numPhases             
            dataPerPhaseAll{n_testPhase} = cell2mat(arrayfun(@(x) mean(x{1,1}(timePhaseLabels == n_testPhase,:)),SessionData, 'UniformOutput', false));
        end 

        for numPhase = 1:numPhases


            if flagLeaveOneOut
                cv = cvpartition(labelsPerSession, 'LeaveOut');
            else
                cv = cvpartition(labelsPerSession, 'KFold', k);
            end 

            for cvRun = 1:cv.NumTestSets %

                trIdx = find(cv.training(cvRun));
                teIdx = find(cv.test(cvRun));

                dataTrainingPhase = dataPerPhaseAll{numPhase}(trIdx,:);
                labelsTrain = labelsPerSession(trIdx); 

                dataTestingPhase = cellfun(@(x) x(teIdx,:), dataPerPhaseAll, 'UniformOutput', false);
                labelsTest = labelsPerSession(teIdx); 

                 for n_testPhase = 1:numPhases

                      [errTrain,errTestCv, ~, predictedTestGeneral] = ...
                    classification.LDA_classification_train_test(dataTrainingPhase, dataTestingPhase{n_testPhase},labelsTrain, labelsTest, ...
                  'flagPCA', flagPCA, 'PCA_variance', PCA_variance, 'classifierType', classifier, 'flagTuning', false); 

                    errTestTmp(cvRun, n_testPhase) = errTestCv;

                 end           
            end 

            errTest(rep,numPhase,:, n_session) = mean(errTestTmp,1);
        end 
    end     
end 

if flagSaveData
     dataFolderPath = fullfile(pwd,'Database',subjectId,'Data','ClassificationData','CrossPhaseClassification');
    save(fullfile(dataFolderPath,[unit_region '_' taskCue  '_Speech_Data']),'errTest');
end 

