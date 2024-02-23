%% compute Internal Strategy results (Figure S5 of manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

% Add current folder and its subfolders to the search path
addpath(genpath(pwd)); 
clc;
clear all;
%close all; 

%%
subjectId = 's2';
unitRegion = 'SMG'; %all, PMV, SMG, AIP, or BA5, or S1X

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd, 'Database', subjectId, 'Data', 'SupplementaryFigure');

cueforSpeech = 'Written_Audio'; % chose one of the four options
%cueforSpeech = 'Written_Written';
%cueforSpeech = 'Audio_Written';
%cueforSpeech = 'Audio_Audio';

cueType = ['Table_Cross_Cue_' cueforSpeech '_perTrials.mat'];

Data = load(fullfile(dataFolderPath, cueType));  
Go_data = Data.Go_data;

sessions = unique(Go_data.session_date);

numSessions = length(sessions);
flagTunedChannels = false; %select tuned channels through krukalwallis test
k = 8;
classifier = 'diaglinear'; 
numRepetition = 1; 
numPhases = 6; 
numChannels = 96;
FontSize = 12;
PCA_variance = 90;
flagPCA = true;
flagLeaveOneOut = true; 
flagSaveData = false;

for sessionIdx = 1:numSessions

    disp(['Classifying session ' sessions{sessionIdx} ]);
    
    if numSessions ~= 1
        idxThisSession = ismember(Go_data.session_date, sessions(sessionIdx));  
    else  
        idxThisSession = ismember(Go_data.session_date, sessions);
    end 
    
    SessionDataAll = Go_data(idxThisSession,:);
    labels = Go_data.GoLabels(idxThisSession);
    timePhaseLabelsAll = Go_data.time_phase_labels(idxThisSession);

    if strcmp('SMG', unitRegion)
        sessionData= SessionDataAll.SMG_Go;
    elseif strcmp('S1X', unitRegion)
        sessionData = SessionDataAll.S1X_Go;
    end
    
    timePhaseLabels =timePhaseLabelsAll{1};
    
    labelsPerSession = labels; 
    

    for rep = 1:numRepetition 
        
        DataPerPhaseAll = cell(numPhases,1);
        
        for testNbr = 1:numPhases
            
            DataPerPhaseAll{testNbr} = cell2mat(arrayfun(@(x) mean(x{1,1}(timePhaseLabels == testNbr,:)),sessionData, 'UniformOutput', false));
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

                dataTraining = DataPerPhaseAll{numPhase}(trIdx,:);
                labelsTrain = labelsPerSession(trIdx); 
                
                dataTesting = cellfun(@(x) x(teIdx,:), DataPerPhaseAll, 'UniformOutput', false);                    
                labelsTest = labelsPerSession(teIdx); 
                                  
                                  
                for testNbr = 1:numPhases
                    [errTrain,errTestCv, ~, predictedTestGeneral] = ...
                    classification.LDA_classification_train_test(dataTraining, dataTesting{testNbr},labelsTrain, labelsTest, ...
                    'flagPCA', flagPCA, 'PCA_variance', PCA_variance, 'classifierType', classifier, 'flagTuning', flagTunedChannels); 
                    
                    errTest(rep,cvRun,numPhase,testNbr, sessionIdx) = errTestCv;

                end           
            end 
        end 
    end  
end 


errTesting = (1- squeeze(mean(errTest,2)))*100;
dataFolderPath = fullfile(pwd, 'Database', subjectId, 'Data', 'ClassificationData', 'CrossPhaseInternalStrategy');     

if flagSaveData
    save(fullfile(dataFolderPath, [cueforSpeech '.mat']),'errTesting' )
end

