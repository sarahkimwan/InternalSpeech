function [tunedCombinedChannels, tunedChannelsPerPhase, tunedChannelsPerBin, sumPhase, sumBin,numTunedChannelsPerCategory, tunedChannelsPerPhasePerCategory,numTunedChannelsPerCategoryBin,pPerPhase, pPerPhaseOrig] = getRegressionTunedChannels(Data,Labels, TimePhaseLabels,varargin)

    % Takes as input a neural cell array [numTrialsx1] with experimental labels [numTrialsx1] and cell array of trial phase labels [numTrialsx1] 
    % returns a bool of tuned units using a linear regression analysis 
    
    %default values
    multipleCompare = true;
    combineTunedChannels = true;
    flagBinperBin = false;
    flagShuffleTest = false; 
    
    % Loading optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'multcompare'
                multipleCompare = varargin{2};
            case 'combinetunedchannels'
                combineTunedChannels = varargin{2};
            case 'binperbintuning'
                flagBinperBin = varargin{2};                          
            case 'flagshuffletest'
                flagShuffleTest = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
        end
          varargin(1:2) = [];
    end
    
    %class names of trials
    categoryNames = utile.image2class_simple(unique(Labels));

    alpha = 0.05;
    numBins = size(Data{1},1);
    
    numConditions = numel(unique(Labels));
    
    [groupCount,~] = groupcounts(Labels);
    
    Condition_names = [{'ITI'}];
    Condition = repmat(Condition_names, [max(groupCount),1]);
    
    ConLabels =utile.image2class_simple(Labels)';
    
    Condition = [Condition; ConLabels];
    
    %Presence vector X
    for i = 1:numConditions
        X(:,i) = ismember(Condition, categoryNames{i});
    end
    
    numChannels = size(Data{1},2);
    numPhases = numel(unique(TimePhaseLabels{1}));
    
    pPerPhase = ones(numPhases,numConditions, numChannels);
    
    pPerPhaseOrig = ones(numPhases, numConditions,numChannels);
    pPerBin = ones(numBins, numConditions, numChannels);
    pFStat = ones(numPhases,numChannels);
    rPerPhase = zeros(numPhases,numChannels); 
    rPerPhaseAdjusted = zeros(numPhases,numChannels); 
    
    pMultcompPerPhase = ones(numPhases,numConditions,numChannels);
    pMultcompPerBin = ones(numBins, numConditions, numChannels);
    ITI_DataAll =  cell2mat(cellfun(@(x,y) mean(x(y== 1,:),1),Data,TimePhaseLabels, 'UniformOutput', false));

  
    for nPhase = 1:numPhases
        
        DataPerPhase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== nPhase,:),1),Data,TimePhaseLabels, 'UniformOutput', false));

    
      for nChannel = 1:size(DataPerPhase,2)
          
          DataPerTrial = DataPerPhase(:,nChannel);
    
          if flagShuffleTest % randomize trials for shuffle distribution
              DataPerTrial = DataPerTrial(randperm(length(DataPerTrial)));
          end
       
          ITI_Data = ones(max(groupCount),1)*mean(mean(ITI_DataAll(:,nChannel)));
          
          %Perform linear regression test for each channel for each phase        
          FR = [ITI_Data; DataPerTrial]; 
          mdl = fitlm(X,FR);
          pPerPhase(nPhase,:,nChannel) = mdl.Coefficients.pValue(2:end); 
          pPerPhaseOrig(nPhase,:,nChannel) =  pPerPhase(nPhase,:,nChannel);
          rPerPhase(nPhase, nChannel) = mdl.Rsquared.Ordinary;
          rPerPhaseAdjusted(nPhase, nChannel) = mdl.Rsquared.Adjusted;
        
      end    
    end 
    
    %perform analysis for each time bin
    if flagBinperBin
    
        for nBin = 1:numBins
            disp(['Bin number ' num2str(nBin)])
            DataPerBin = cell2mat(arrayfun(@(x,y) Data{x,1}(nBin,:),1:size(Data,1), 'UniformOutput', false)');
    
            for nChannel = 1:numChannels
    
              DataPerBinTrial = DataPerBin(:,nChannel);
              DataPerBinTrialOrdered = DataPerBinTrial;
              ITI_Data = ones(max(groupCount),1)*mean(mean(ITI_DataAll(:,nChannel)));
    
              %perform linear regression for each channel for each phase for each of the 5 grasps 
              FR = [ITI_Data; DataPerBinTrialOrdered];  
              mdl = fitlm(X,FR);
              p_val_bin = mdl.Coefficients.pValue(2:end);
              
              if multipleCompare
                  [~,~,p_val_bin] = utile.MultipleComparisonsCorrection(p_val_bin,'method', 'fdr');
              end 
          
              pPerBin(nBin,:,nChannel) = p_val_bin; 
    
            end
        end 
    end 
    
    
    if multipleCompare
        disp('Performing multiple comparison')
    
        for nChannel = 1:numChannels
            [~,~,pMultcompPerPhase(:,:,nChannel)] = utile.MultipleComparisonsCorrection(pPerPhase(:,:,nChannel),'method', 'fdr'); 
        end 
        %replace old p values by multcompare p values
        pPerPhase = pMultcompPerPhase; 
    end 
    
    
    tunedChannelsPerPhasePerCategory = pPerPhase < alpha;
    numTunedChannelsPerCategory = sum(tunedChannelsPerPhasePerCategory,3);
    numTunedUnitsPhase = cell2mat(arrayfun(@(x) nnz(squeeze(sum(squeeze(tunedChannelsPerPhasePerCategory(x,:,:))))), 1:size(tunedChannelsPerPhasePerCategory,1), 'UniformOutput', false));
     
    tunedChannelsPerPhase = logical(squeeze(sum(tunedChannelsPerPhasePerCategory,2)))';
    sumPhase = sum(tunedChannelsPerPhase);
    FStatSumPhase = sum(pFStat < 0.05,2);
    tunedChannelsPerBin = pPerBin < 0.05;
    numTunedChannelsPerCategoryBin = sum(tunedChannelsPerBin,3);
    
    sumBin = sum(logical(squeeze(sum(tunedChannelsPerBin,2))),2); 
    
    if combineTunedChannels
        tunedCombinedChannels = logical(sum(tunedChannelsPerPhase,2));    
    end

end

