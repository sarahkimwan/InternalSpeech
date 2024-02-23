%% plot offline decoding results (Figure 5A,C and S4 of manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

% Add current folder and its subfolders to the search path
addpath(genpath(pwd)); 

clc;
clear all;

% Define the subject ID
subjectId = 's2'; % for subject s2
%subjectId = 's3'; % for subject s3
unitRegion = 'SMG'; % SMG or S1X


%% Parameters

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd, 'Database', subjectId, 'Data', 'ClassificationData', 'ClassificationPerTrial');
dataName = fullfile(dataFolderPath, [unitRegion '_All_Speech_Data']);
shuffleName = fullfile(dataFolderPath, [unitRegion '_All_Speech_Shuffle']);

% Load data
data = load(dataName);
shuffle = load(shuffleName);

% Define phase names
phaseNames = {'ITI', 'Cue', 'D1', 'Imagined', 'D2', 'Speech'};
numPhases = length(phaseNames);

% Calculate mean accuracy
keepMeanAcc = (1 - squeeze(nanmean(data.errTest, 2))) * 100;
numSessions = size(keepMeanAcc, 2);

% Identify and exclude empty sessions
emptySessions = find(isnan(keepMeanAcc(1, :)));
sessionsToKeep = setdiff(1:numSessions, emptySessions);

% Process shuffle data
shuffleData = (1 - squeeze(nanmean(shuffle.errTest, 2))) * 100;
keepMeanAcc = keepMeanAcc(:, sessionsToKeep);
shuffleData = shuffleData(:, :, sessionsToKeep);

% Reshape shuffle data for analysis
shuffleData = permute(shuffleData, [1, 3, 2]);
shuffleDataR = reshape(shuffleData, [], numPhases);

% Calculate shuffle mean and percentile values
shuffleMean = squeeze(nanmean(shuffleData, 1));
prcVal = prctile(shuffleDataR, [97.5, 99.5]);

% Cohen's d for shuffle distribution 
for n_phase = 1:numPhases
    cohens_d =utile.cohens_d(keepMeanAcc(n_phase, :), shuffleDataR(n_phase, :), 'independent');
    disp([phaseNames{n_phase} ' - Shuffle: cohens_d ' num2str(cohens_d,3)]);
end

% Plot classification results
fig = figure();
colorModel = [62, 164, 144] / 255;

% Calculate error for confidence interval
errCi = utile.calculateErrCi(keepMeanAcc);

disp(['Confidence interval of the mean ' num2str(errCi(1,:),3)])

plot(squeeze(keepMeanAcc), 'Marker', '.', 'MarkerSize', 12, 'LineStyle', 'none', 'Color', [0, 0, 0]);
hold on
plot(nanmean(shuffleMean, 2), 'Marker', '.', 'MarkerSize', 15, 'LineStyle', '--', 'Color', 'red')
hold on

errorbar(1:size(keepMeanAcc, 1), nanmean(keepMeanAcc, 2), errCi(1, :), errCi(2, :), '-s', 'LineStyle', '--', 'MarkerSize', 10, 'MarkerEdgeColor', colorModel, 'MarkerFaceColor', colorModel, 'Color', colorModel) 

% Set x-axis ticks and labels
xticks(1:numPhases);
xticklabels(phaseNames)
xtickangle(-45)

ylabel('Classification accuracy [%]');
title('Classification combining data');
fontSize = 15;

% Add mean values as text on the plot
meanAll = mean(keepMeanAcc, 2);
for nPhase = 1:numPhases
    text(nPhase + 0.2, meanAll(nPhase), [num2str(meanAll(nPhase), 2) '%'], 'FontSize', fontSize)
end 

% Determine maximum accuracy for setting plot limits
maxAcc = max(keepMeanAcc') + 4;
valAdd = 3;
xlim([0.5, 6.5])

%%
% Statistical comparison between classification accuracy in different phases

font_size = 12;
pVal_all = [];
for nPhaseToTest = [1,3,5]
    [~, pVal] = ttest(keepMeanAcc(nPhaseToTest, :), keepMeanAcc(nPhaseToTest+1, :));
    cohens_d =utile.cohens_d(keepMeanAcc(nPhaseToTest, :), keepMeanAcc(nPhaseToTest+1, :), 'paired');

    pVal_all = vertcat(pVal_all,  pVal);
    if pVal < 0.05
        text(nPhaseToTest + 0.5,maxAcc(nPhaseToTest+1) + valAdd +8, return_sign(pVal), 'Color', 'black', 'LineWidth', 2,'HorizontalAlignment','center','FontSize', font_size);
        disp([phaseNames{nPhaseToTest} ' vs. ' phaseNames{nPhaseToTest+1}])
        disp(['cohens_d ' num2str(cohens_d,3)]);
        utile.sigline([nPhaseToTest-0.1625,nPhaseToTest+1+0.1625],'', [], maxAcc(nPhaseToTest+1) + valAdd);
    end
end

%%
% Highlight phases with significanlty higher classificationa accuracy
% comped to shuffle distribution 
sigVal = prcVal < mean(keepMeanAcc, 2)';

for n_phase = 1:numPhases
    
     if sigVal(2,n_phase)
        text(n_phase,maxAcc(n_phase) , '**', 'Color', colorModel, 'LineWidth', 10, 'FontSize', font_size,'HorizontalAlignment','center'); 
     elseif sigVal(1,n_phase)
         text(n_phase,maxAcc(n_phase) , '*', 'Color', colorModel, 'LineWidth', 10, 'FontSize', font_size,'HorizontalAlignment','center'); 
     end
    
  
end 
set(gca, 'FontSize', fontSize)

if strcmp(unitRegion, 'S1X') || strcmp(subjectId, 's3')
    ylim([0, 70]) 
else
    ylim([0, 120]) 
end
yticks([0:20:100])

% Plot classification matrix
realLabels = squeeze(data.keepTestingLabels);
predLabels = squeeze(data.keepPredictedLabels);
errMat = figure('units', 'normalized', 'outerposition', [0, 0, 0.4, 0.8]);

for phaseIndex = 1:numPhases
    subplot(3, 2, phaseIndex)
    % Combine labels from all session days
    realLabelsComb = reshape(realLabels(:, phaseIndex, :), [], 1);
    predLabelsComb = reshape(predLabels(:, phaseIndex, :), [], 1);

    % Remove NaN values
    realLabelsComb = realLabelsComb(~isnan(realLabelsComb));
    predLabelsComb = predLabelsComb(~isnan(predLabelsComb));

    % Create confusion matrix
    [confusionMat, groupNames] = confusionmat(realLabelsComb, predLabelsComb);
    colormap(jet)
    imagesc(confusionMat ./ sum(confusionMat, 2) * 100)

    % find condition names
    condToTest = utile.image2class_simple(unique(realLabelsComb));
     %plot axis
    if  ismember(n_phase,[5,6])
        xticks([1:8]);
        xticklabels(condToTest)
        xtickangle(-45)
        xlabel('Predicted Class')

    else
        set(gca,'XTick',[])
    end 
    
    %
    if ismember(n_phase,[1,3,5])
        yticks([1:8])
        yticklabels(condToTest)
        ylabel('True Class')
    else
        set(gca,'YTick',[])
    end 
    caxis([0 100])

    title([' Phase : ' phaseNames{n_phase}] )
    set(gca, 'fontSize', 12)
end



%% Helper Functions

function sigSign = return_sign(p)

    if p < 0.001
        sigSign = '***';
    elseif p < 0.01
        sigSign = '**';
    elseif p <= 0.05
        sigSign = '*';
    elseif p > 0.05
        sigSign = '';
    end 

end 

