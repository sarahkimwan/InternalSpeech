function [errTrain,errTest,predictedTrain, predictedTest] = LDA_classification_train_test(dataTrain, dataTest, labelsTrain, labelsTest, varargin)
% Perform LDA classification based on a input training and testing dataset

    %should be default but can be changed
    flagPCA = true;
    PCA_variance = 90;
    classifierType = 'diaglinear';
    flagTuning = false;
    
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
            otherwise
                error(['Unexpected option: ' varargin{1}])
        end
          varargin(1:2) = [];
    end
    % select features calculating tuned units 
    if flagTuning
        numChannels = size(dataTrain,2);
        p_val = nan(1,numChannels);
        for n_Channel = 1:numChannels
            dataPerTrial = dataTrain(:,n_Channel);
            [p_val(n_Channel), ~, ~] = kruskalwallis(dataPerTrial,labelsTrain,'off');
        end
        
        %keep channels that are tuned (alpha = 0.05)
        tunedChannels = p_val < 0.05; 
    
        %if no enough tuned channels available - keep all
        if nnz(tunedChannels) == 0
            tunedChannels = 1:numChannels;
            warning('no tuned channels available - taking all channels ')
        end 
        
        dataTrain = dataTrain(:,tunedChannels);
        dataTest = dataTest(:,tunedChannels); 
    end  
                
    % select features by performing PCA
    if flagPCA
        %calculate pca coefficient on training data
        [coeff, ~, ~,~, explained] = pca(dataTrain); 
        %calculate explained variance
        variance = cumsum(explained); 
        %find number of PCs explaining given threshold
        PCAidx = find(variance > PCA_variance,1); 
        %project data on coefficients
        PCAdataTrain = dataTrain*coeff;
        PCAdataTest = dataTest*coeff;
        dataTrain = PCAdataTrain(:,1:PCAidx);
        dataTest = PCAdataTest(:,1:PCAidx);
    end 
    
    %Train model
    model = fitcdiscr(dataTrain, labelsTrain, 'DiscrimType', classifierType);
    
    %Predict train and test labels
    predictedTrain = predict(model, dataTrain);
    predictedTest = predict(model, dataTest);
    
    %calculate error
    errTrain = 1-(nnz(labelsTrain == predictedTrain)/numel(labelsTrain));
    errTest = 1-(nnz(labelsTest == predictedTest)/numel(labelsTest));

end





