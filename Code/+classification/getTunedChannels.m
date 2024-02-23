function [TunedCombinedChannels, TunedChannelsPerPhase, TunedChannelsPerBin, sumPhase, sumBin] = getTunedChannels(Data,Labels, TimePhaseLabels,varargin)

    % Takes as input a neural cell array [numTrialsx1] with experimental labels [numTrialsx1] and cell array of trial phase labels [numTrialsx1] 
    % returns a bool of tuned units using the kruskalwallis test
    
    %default values
    flagMultipleCompare = true;
    flagBinperBin = false;
    flagCombineTunedChannels = true;
    flagRemoveITItuning = false;
    flagShuffleTest = false;
    
    
    % Loading optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'multcompare'
                flagMultipleCompare = varargin{2};
            case 'binperbintuning'
                flagBinperBin = varargin{2};              
            case 'combinetunedchannels'
                flagCombineTunedChannels = varargin{2};             
            case 'removeitituning'
                flagRemoveITItuning = varargin{2};
            case 'flagshuffletest'
                flagShuffleTest = varargin{2};           
            otherwise
                error(['Unexpected option: ' varargin{1}])
        end
          varargin(1:2) = [];
    end
    
    
    numChannels = size(Data{1},2);
    uniqueTmp = unique(TimePhaseLabels{1});
    numPhases = length(uniqueTmp);
    numBins = size(Data{1},1);
    numUnits = size(Data{1},2);
    pPerPhase = ones(numChannels, numPhases);
    pPerBin = ones(numChannels, numBins);
    
    pMultcompPerPhase = ones(numChannels, numPhases);
    pMultcompPerBin = ones(numChannels, numBins);
    
    for n_phase = 1:numPhases
      
        DataPerPhase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== uniqueTmp(n_phase),:),1),Data,TimePhaseLabels, 'UniformOutput', false));
            
     if flagShuffleTest
         Labels = Labels(randperm(length(Labels)));
     end 
        
      for n_channel = 1:numChannels
          
          DataPerTrial = DataPerPhase(:,n_channel);
          %perform kruskal wallis test 
          [pPerPhase(n_channel,n_phase), ~, ~] = kruskalwallis(DataPerTrial,Labels,'off');
    
      end
      
           
    end 
    
    if flagBinperBin
    
        for n_bin = 1:numBins
            disp([ 'Bin nbr ' num2str(n_bin)]);
            DataPerBin = cell2mat(arrayfun(@(x,y) Data{x,1}(n_bin,:),1:size(Data,1), 'UniformOutput', false)');
    
            for n_channel = 1:numChannels
    
              DataPerBinTrial = DataPerBin(:,n_channel);
              %perform kruskal wallis test 
              [p_val_bin, ~, ~] = kruskalwallis(DataPerBinTrial,Labels,'off');
              
              pPerBin(n_channel,n_bin) = p_val_bin;
    
            end
    
        end 
    
    end 
    
    %adjust p-values for multiple comparisons
    
    if flagMultipleCompare
    
        for n_channel = 1:numChannels
            [~,~,pMultcompPerPhase(n_channel,:)] = utile.MultipleComparisonsCorrection(pPerPhase(n_channel,:),'method', 'fdr'); 
        end 
    
        %replace old p values by multcompare p values
        pPerPhase = pMultcompPerPhase; 
    
    end 
    
    TunedChannelsPerPhase = pPerPhase < 0.05;
    sumPhase = sum(TunedChannelsPerPhase);
    
    TunedChannelsPerBin = pPerBin < 0.05;
    sumBin = sum(TunedChannelsPerBin); 
    % Combine the units that are tuned in either of the phases  - can be used
    % for classification
    if flagCombineTunedChannels
        if flagRemoveITItuning
            TunedCombinedChannels = logical(sum(pPerPhase(:,2:end) < 0.05,2));
        else
            TunedCombinedChannels = logical(sum(pPerPhase < 0.05,2));
        end 
        
    end 

end

