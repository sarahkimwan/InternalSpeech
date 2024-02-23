function [errTrain,errTest, labelsTestAll,predictedTestAll] = LDA_classification_rep(data,labels, varargin)

    %perform cross validation LDA classification based on 'data' [trials  x features]
    % and 'labels' [trials x 1]
    
    %default values
    flagPCA = true;
    PCA_variance = 90;
    classifierType = 'diaglinear';
    flagTuning = false;
    k = 8; 
    flagLeaveOneOut = true; 
    flagRandomPerm = false; 
    numRep = 1;
    flagErrorMatrix = false;
    
    
    % Loading optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'flagpca'
                flagPCA = varargin{2};
            case 'pca_variance'
                PCA_variance = varargin{2};
            case 'classifiertype'
                classifierType = varargin{2};              
            case 'flagtuning'
                flagTuning = varargin{2};
            case 'k'
                k = varargin{2};               
            case 'flagleaveoneout'
                flagLeaveOneOut = varargin{2};
            case 'flagrandomperm'
                flagRandomPerm = varargin{2};           
            case 'numrep'
                numRep = varargin{2};
            case 'flagerrormatrix'
                flagErrorMatrix = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
        end
          varargin(1:2) = [];
    end
    
    
    if flagTuning && flagPCA
        error('Do you want to do both tuning and PCA?')
        %if on purpose remove 
    end 
    
    
    if size(data,1) ~= size(labels,1)
        error('Wrong label size')
    end
    
    if flagLeaveOneOut
        numFolds = size(data,1);
    else
        numFolds = k; 
    end
    
    errTrain = nan*ones(numRep,numFolds);
    errTest = nan*ones(numRep,numFolds);
    labelsTestAll = nan*ones(numRep,size(data,1));
    predictedTestAll = nan*ones(numRep,size(data,1));
    
    %iterate over number of repetitions
    for rep = 1:numRep
        
        if numFolds == size(data,1)
            cv = cvpartition(size(data,1),'LeaveOut');
        else
            cv = cvpartition(size(data,1), 'KFold', k);
        end 
        
        %randomize trial labels
        if flagRandomPerm
            labels = labels(randperm(length(labels)));
        end 
        
        labelsTestTmp = [];
        predictedTestTmp = [];
        
        for cvRun = 1:cv.NumTestSets %
    
            trIdx = find(cv.training(cvRun));
            teIdx = find(cv.test(cvRun));
            %training set
            dataTrain = data(trIdx,:);
            labelsTrain = labels(trIdx); 
            %testing set
            dataTest = data(teIdx,:);        
            labelsTest = labels(teIdx);
            labelsTestTmp = [labelsTestTmp; labelsTest];
    
            if flagTuning
                numChannels = size(dataTrain,2);
                p_val = nan(1,numChannels);
                for ChannelNbr = 1:numChannels
    
                    DataPerTrial = dataTrain(:,ChannelNbr);
                     [p_val(ChannelNbr), ~, ~] = kruskalwallis(DataPerTrial,labelsTrain,'off');
    
                end
                
                tunedChannels = p_val < 0.05;
                
                if nnz(tunedChannels) == 0
                    tunedChannels = 1:numChannels;
                    warning('no tuned channels available - taking all channels ')
                end 
                
                %disp('# tuned channels ')
               % nnz(tunedChannels)
                dataTrain = dataTrain(:,tunedChannels);
                dataTest = dataTest(:,tunedChannels); 
                           
            end 
            
            %perform PCA
            if flagPCA
                [coeff, ~, ~,~, explained] = pca(dataTrain); 
                
                variance = cumsum(explained); 
                %select number of PCs based on percentage of variance determined by PCA_variance
                PCA_idx = find(variance > PCA_variance,1); 
                PCA_DataTrain = dataTrain*coeff;
                PCA_DataTest = dataTest*coeff;
                dataTrain = PCA_DataTrain(:,1:PCA_idx);
                dataTest = PCA_DataTest(:,1:PCA_idx);
            end 
            
            model = fitcdiscr(dataTrain, labelsTrain, 'DiscrimType', classifierType);
    
            predictedTrain = predict(model, dataTrain);
            predictedTest = predict(model, dataTest);
            
            predictedTestTmp = [predictedTestTmp;predictedTest];
            
            errTrain(rep,cvRun) = 1-(nnz(labelsTrain == predictedTrain)/numel(labelsTrain));
            errTest(rep,cvRun) = 1-(nnz(labelsTest == predictedTest)/numel(labelsTest));
                
            clear labelsTest
            clear predictedTest
        end
        
        labelsTestAll(rep,:) = labelsTestTmp;
        predictedTestAll(rep,:) = predictedTestTmp;
    
    end 
    
    %plot error matrix
    if flagErrorMatrix
        figure();
        labelsTest = labelsTestAll;
        
        predictedTestAll = predictedTestAll;
        
        [cm,gn] = confusionmat(labelsTest, predictedTestAll);
        class_names = utile.image2class_simple(gn);
    
        confusionchart(cm,class_names)
    end


end

