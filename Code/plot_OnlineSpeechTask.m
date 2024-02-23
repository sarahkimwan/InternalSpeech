%% plot online decoding results (Figure 5B and 5C of manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); %add folder to search path 
clc
clear all
%close all

%% parameters

subject_id = 's2'; %define subject_id
%subject_id = 's3'; %define subject_id
unit_region = 'SMG'; 

dataFolderPath = fullfile(pwd,'Database', subject_id,'Data','ClassificationData','Online');

%load data
Data = load(fullfile(dataFolderPath,'Data'));
Shuffle = load(fullfile(dataFolderPath,'Shuffle_1000'));
trainingDatapointsPerRun = Data.trainingDatapointsPerRun;
decodingAccuraciesAll = Data.decodingAccuraciesAll;
decodingAccuraciesShuffle = Shuffle.decodingAccuraciesShuffle;

if strcmp(subject_id, 's3')
    %For s3: used first two datablocks for training
    trainingDatapointsPerRun = cellfun(@(x) x(2:end), trainingDatapointsPerRun, 'UniformOutput', false);
    decodingAccuraciesShuffle = cellfun(@(x) x(:,2:end), decodingAccuraciesShuffle, 'UniformOutput', false);
    signMax = 40;
    
else
    signMax = 100;
end 

% online decoding results
onlineDecoding = cell2mat(decodingAccuraciesAll); 
% number of training trials used to train the model for each run
nTrainingDatapointsPerRun = cell2mat(trainingDatapointsPerRun); 
% shuffle results per run
shufflePerRun = cell2mat(cellfun(@(x) mean(cell2mat(x)), decodingAccuraciesShuffle, 'UniformOutput', false)); 
% shuffle results concatenated over all runs
shuffleAll = cell2mat(cellfun(@(x) cell2mat(x), decodingAccuraciesShuffle, 'UniformOutput', false)); 
% sort results based on number of training trials
[trainingSorted,b] = sort(nTrainingDatapointsPerRun);
decodingSorted = onlineDecoding(b)*100; 


shuffleSorted = shufflePerRun(b); 
shuffleAllSorted = shuffleAll(:,b);
categories = [];

for nRun = 1:length(trainingSorted)
    %categorize runs according to more or less than 16 trials for training
    if trainingSorted(nRun) < 16
        categories(nRun) = 1;
    elseif trainingSorted(nRun) >= 16
        categories(nRun) = 2;
    end 
end

%sort shuffle according to category
shuffleLow = reshape(shuffleAllSorted(:,categories ==1),[],1);
shuffleHigh = reshape(shuffleAllSorted(:,categories ==2),[],1);

%define percentiles to define significance 
prcLow = prctile(shuffleLow, [97.5, 99.5, 99.95]);
prcHigh = prctile(shuffleHigh, [97.5, 99.5, 99.95]);

lowDec = decodingSorted(categories ==1);
highDec = decodingSorted(categories ==2);
%calclate confidence intervals 
errCi(:,1) = utile.calculateErrCi(lowDec);
errCi(:,2) = utile.calculateErrCi(highDec);

%% plot figure 5B of manuscript
fig = figure(); 
plot(categories, decodingSorted, 'Marker', '.', 'MarkerSize', 12, 'LineStyle','none', 'Color', [0,0,0]);
hold on
plot(categories, shuffleSorted, 'Marker', '.', 'MarkerSize', 12, 'LineStyle','--', 'Color', [1,0,0]);

ylim([0 110])
xlim([0.5 2.5])
xticks([1:2]);
if strcmp(subject_id, 's2')
    xticklabels({'8 - 14', '16 - 20'})
else
    xticklabels({'8 - 14', '16 - 32'})
end 
ylabel('Online decoding accuracy')
xlabel('Number of training trials per word')
hold on

font_size = 12;
ColorModel = [0 0 1];
err = errorbar((1:2), [mean(lowDec), mean(highDec)], errCi(1,:), errCi(2,:) ,'-s', 'LineStyle', 'none','MarkerSize',12,...
             'MarkerEdgeColor',ColorModel,'MarkerFaceColor',ColorModel , 'Color', ColorModel);
      
meanAll = [mean(lowDec), mean(highDec)];

for nCategory = 1:2
    text(nCategory+0.1, meanAll(nCategory),[num2str(meanAll(nCategory),2) '%'],'FontSize', font_size)
end 

%plot percentiles 
sigLow =  mean(lowDec) > prcLow;
sigHigh = mean(highDec) > prcHigh;

if  nnz(sigLow)
    text(1,signMax, returnSign(sigLow), 'Color', 'blue', 'LineWidth', 2,'HorizontalAlignment','center','FontSize', font_size);
end

if  nnz(sigHigh)
    text(2,signMax, returnSign(sigHigh), 'Color', 'blue', 'LineWidth', 2,'HorizontalAlignment','center','FontSize', font_size);
end

%calcualte significance between low and high number of training trials for
%decoding accuracies 
[h,p] = ttest2(decodingSorted(categories ==1), decodingSorted(categories == 2));

%plot values
if p < 0.05
    text(1.5,100, ['*'], 'Color', 'black', 'LineWidth', 2,'HorizontalAlignment','center','FontSize', font_size);
    text(1.5,90, [ 'p = ' num2str(p,2) ], 'Color', 'black', 'LineWidth', 2,'HorizontalAlignment','center','FontSize', font_size);
    utile.sigline([1,2],'', [], 92);
end
set(gca, 'FontSize', font_size)



%% plot confusion matrix (Figure 5C)

wordLabels ={'Battlefield', 'Cowboy', 'Python', 'Spoon', 'Swimming', 'Telephone', 'Bindip', 'Nifzig'};

tableData = Data.table2Name;

labelNames = tableData.LabelNames;
Id = tableData.Id;

adaptOrder = cell2mat(cellfun(@(x) find(ismember(labelNames,x)), wordLabels, 'UniformOutput', false));
labelsCued = Data.LabelsCued; 
labelsPredicted = Data.LabelsPredicted;

if strcmp(subject_id, 's2')
    LabelsA1 = cell2mat(labelsCued(2:end,1));
    LabelsA2 = cell2mat(labelsCued(2:end,2));
    LabelsA3 = cell2mat(labelsCued(2:end,3));
else
    LabelsA1 = cell2mat(labelsCued(3:end,1));
    LabelsA2 = cell2mat(labelsCued(3:end,2));
    LabelsA3 = cell2mat(labelsCued(3:end,3));
end 


labelsP1 = cell2mat(labelsPredicted(2:end,1));
labelsP2 = cell2mat(labelsPredicted(2:end,2));
labelsP3 = cell2mat(labelsPredicted(2:end,3));

labelsCuedConcat = [LabelsA1;LabelsA2;LabelsA3];
labelsPredictedConcat = [labelsP1;labelsP2;labelsP3];

confusionMat = figure(); 
colormap(jet)

[C,g] = confusionmat(labelsCuedConcat, labelsPredictedConcat, 'Order',Id(adaptOrder));

number_classes = length(wordLabels);
imagesc(C./sum(C,2)*100)
xticks([1:number_classes])
xticklabels(wordLabels)
xtickangle(45)
yticks([1:number_classes])
yticklabels(wordLabels)
xlabel('Predicted Class')
ylabel('True Class')
caxis([0 100])

%% function

function sigSign = returnSign(p)

    if p(3)
        sigSign = '***';
    elseif p(2)
        sigSign = '**';
    elseif p(1)
        sigSign = '*';
    else
        sigSign = '';
    end 

end 
