%% extract and compute data for online classification plots

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway
% Add current folder and its subfolders to the search path
addpath(genpath(pwd)); 

clc
clear all
%close all
%%

subjectId = 's2'; %or s3 participant
%subjectId = 's3';

flagSavePlot = false;
flagShuffle = false; %compute shuffle distribution 

if strcmp(subjectId, 's2')

taskfiles_all{1} = {'20220817-125234-125244-SpeechTrainingInternal.mat',...
    '20220817-130208-130644-SpeechTrainingInternal.mat',...
    '20220817-131355-131634-SpeechTrainingInternal.mat',...
    '20220817-132801-132822-SpeechTrainingInternal.mat',...
  };

% session day online 20221107 

taskfiles_all{2} = {'20221107-125652-125704-SpeechTrainingInternal.mat',...
    '20221107-130618-131258-SpeechTrainingInternal.mat',...
    '20221107-132245-132532-SpeechTrainingInternal.mat',...
    '20221107-133115-133152-SpeechTrainingInternal.mat',...
  };

% session day online 20221212 
taskfiles_all{3} = ...
   {'20221212-115426-115443-SpeechTrainingInternal.mat',...
    '20221212-120852-120901-SpeechTrainingInternal.mat',...
    '20221212-121219-121405-SpeechTrainingInternal.mat',...
    '20221212-121901-121959-SpeechTrainingInternal.mat',...
    '20221212-122310-122358-SpeechTrainingInternal.mat',...
    '20221212-122649-122720-SpeechTrainingInternal.mat',...
    '20221212-123030-123119-SpeechTrainingInternal.mat',...
    '20221212-123433-123601-SpeechTrainingInternal.mat'
  };

elseif strcmp(subjectId, 's3')
    
    taskfiles_all{1} = ...
   {'20230810-102731-102739-SpeechTrainingInternal.mat',... % training
    '20230810-103657-104513-SpeechTrainingInternal.mat',... % training
    '20230810-110031-110257-SpeechTrainingInternal.mat',... % testing 
    '20230810-111334-111514-SpeechTrainingInternal.mat',... % testing
    '20230810-112832-112856-SpeechTrainingInternal.mat'     % testing
  };

 taskfiles_all{2} = ...
   {'20230824-100727-101137-SpeechTrainingInternal.mat',... % training
    '20230824-102205-102509-SpeechTrainingInternal.mat',... % training
    '20230824-103523-103619-SpeechTrainingInternal.mat',... % testing 
    '20230824-103523-104741-SpeechTrainingInternal.mat',... % testing
  };

 taskfiles_all{3} = ...
   {'20230831-101314-101458-SpeechTrainingInternal.mat',... % training
    '20230831-104443-104721-SpeechTrainingInternal.mat',... % training
    '20230831-110142-110157-SpeechTrainingInternal.mat',... % testing 
    '20230831-111225-111336-SpeechTrainingInternal.mat',... % testing
  };
 
end


meanFRNeuralDataAll = {};
labelsAllTmp = {};
labelsPredicted = {};
errTest = [];
neuralDataOnlineAll = {};

training_datapoints_all = {};
decoding_accuracies_all = {};
decoding_accuracies_could_have_all = {};

for n_day = 1:length(taskfiles_all)
        classification_accuracy_within_task ={};

        taskfiles = taskfiles_all{n_day};
        n_blocks = length(taskfiles); 

    for n_taskfile= 1:length(taskfiles)

        task = load(fullfile(pwd,'Database',subjectId,'Data','OnlineData', taskfiles{n_taskfile}));
        task = task.task;

        neuralDataOnline = task.neuralDataSpeech;

        classification_accuracy_within_task{n_taskfile}= mean([task.trialdata.tr_response] == [task.trialdata.tr_target]);

        responses = vertcat(task.trialdata.ex_success);
        labelsAllTmp{n_taskfile} = [task.trialparams.Image_code]';     

        succ = mean(vertcat(task.trialdata.ex_success));
        if n_taskfile~= 1 && strcmp(subjectId, 's2')
            labelsPredicted{n_taskfile,n_day} = vertcat(task.trialdata.tr_response);
            labelsCued{n_taskfile,n_day} = labelsAllTmp{n_taskfile};

        end 
        
        if n_taskfile~= 1 && n_taskfile~= 2 && strcmp(subjectId, 's3')
            labelsPredicted{n_taskfile,n_day} = vertcat(task.trialdata.tr_response);
            labelsCued{n_taskfile,n_day} = labelsAllTmp{n_taskfile};

        end 

        neuralDataOnlineTmp = cellfun(@(x) x(:,:), neuralDataOnline, 'UniformOutput', false);
        meanFRNeuralData = cell2mat(cellfun(@(x) nanmean(x,2), neuralDataOnline, 'UniformOutput', false));


        neuralDataOnlineAll{n_taskfile} = neuralDataOnlineTmp;
        meanFRNeuralDataAll{n_taskfile} = meanFRNeuralData;

    end 
    
    labelNames = {task.trialparams.image};
    labelsCode = [task.trialparams.Image_code];

    [Id,Code_idx] = unique(labelsCode);
    labelNames = labelNames(Code_idx)';
    
    table2Name = table(labelNames, Id');

    class_acc = cell2mat(classification_accuracy_within_task);
    
    %saved decoding accuracies from decoder 
    if strcmp(subjectId, 's2')
        decoding_accuracies_all{n_day} = class_acc(2:end);

    elseif strcmp(subjectId, 's3')
        decoding_accuracies_all{n_day} = class_acc(3:end);
    end 

    %% produce shuffle distribution

    % compute data for training and testing
    neuralDataTimepointsAll = [];
    for n_taskfile = 1:n_blocks
        neuralDataOnline = neuralDataOnlineAll{n_taskfile};
        neuralDataOnlineTmp = cellfun(@(x) x(:,1:end), neuralDataOnline, 'UniformOutput', false);

        neuralDataTimepointsAll = [neuralDataTimepointsAll,neuralDataOnlineTmp];
        meanFRNeuralData = cell2mat(cellfun(@(x) mean(x,2), neuralDataOnlineTmp, 'UniformOutput', false));
        meanFRNeuralDataAll{n_taskfile} = meanFRNeuralData;

    end 


    errTestAll = {};
    PCA_val = 95;
    if flagShuffle 
        num_rep = 1000; 
    else
        num_rep = 1;
    end 

    predictedTestAll = [];
    labelsTestAll = [];
    training_datapoints = [];

    for n_shuffle = 1:num_rep
        
        for n_taskfile = 1:(n_blocks-1)

            dataTrain = cell2mat(meanFRNeuralDataAll(1:n_taskfile))';
            labelsTrain = cell2mat(labelsAllTmp(1:n_taskfile)');

            DataAll = [dataTrain,labelsTrain];

            if flagShuffle
                labelsTrain = labelsTrain(randperm(length(labelsTrain)));
            end

            dataTest = meanFRNeuralDataAll{n_taskfile+1}';
            labelsTest = cell2mat(labelsAllTmp(n_taskfile+1))';
            
            if n_taskfile == 1
                training_datapoints(n_taskfile) = unique(histc(labelsTrain, unique(labelsTrain)));
            else
                training_datapoints(n_taskfile) = size(dataTest,1)/8 + training_datapoints(n_taskfile-1);
            end

            %PCA 
            [coeff, ~, ~,~, explained] = pca(dataTrain); 
            variance = cumsum(explained); 
            idx90 = find(variance > PCA_val,1); 
            PCA_DataTrain = dataTrain*coeff;
            PCA_DataTest = dataTest*coeff;
            dataTrain = PCA_DataTrain(:,1:idx90);
            dataTest = PCA_DataTest(:,1:idx90);
            
            model = fitcdiscr(dataTrain, labelsTrain, 'DiscrimType', 'diaglinear');
            predictedTrain = predict(model, dataTrain);
            predictedTest = predict(model, dataTest);
            predictedTestAll = [predictedTest; predictedTestAll];
            labelsTestAll = [labelsTest,labelsTestAll];
            errTrain = 1-nnz(labelsTrain == predictedTrain)/numel(labelsTrain);
            errTest = 1-nnz(labelsTest' == predictedTest)/numel(labelsTest);


            errTestAll{n_shuffle,n_taskfile} = (1-errTest)*100; 
           
        end 
    end 
    
    training_datapoints_all{n_day} = training_datapoints;

    if flagShuffle
        %save testing accuracy - for onlined decoding, keep decoding
        %accuracy saved during task
        decoding_accuracies_all{n_day} = errTestAll;
    end
end
% 
if flagSavePlot
    if flagShuffle
        save(['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subjectId '\Data\ClassificationData\Online\Shuffle_' num2str(num_rep)],'training_datapoints_all','decoding_accuracies_all');
    else
        save(['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subjectId '\Data\ClassificationData\Online\Data'],'training_datapoints_all','decoding_accuracies_all', 'labelsCued', 'labelsPredicted', 'table2Name');
    end 
end