%% Plot Kruskal-Wallis and regression tuning plots per phase (Figure 3B,E and S3C,D of manuscript)

%% Important: Run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); % Add folder to search path 
clc
clear all
% close all

%% Parameters

subjectId = 's2'; % s2 or s3
unitRegion = 'SMG'; % SMG or S1X
taskCue = {'Written'}; % Choose cue type: Written or Auditory 
%taskCue = {'Auditory'}; % Choose cue type: Written or Auditory 

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd, 'Database', subjectId, 'Data', 'TuningData');

if strcmp(subjectId, 's2')
    cueIdx = find(ismember({'Auditory', 'Written'}, taskCue));
elseif strcmp(subjectId, 's3') && strcmp(taskCue, 'Written')
    cueIdx = 1; 
elseif strcmp(subjectId, 's3') && ~strcmp(taskCue, 'Written')
    error('No Auditory dataset for subject s3')
end

regSub = {'Regression', 'KruskalWallis'};
numMethods = numel(regSub);
fontSize = 12; 

%% Plot all sessions with 95% confidence interval

h = [];

plottingColorsAll = utile.get_color_rgb_codes({'Speech', 'Written'});

if strcmp(taskCue, 'Auditory')
    plottingColors = plottingColorsAll{1};
elseif strcmp(taskCue, 'Written')
    plottingColors = plottingColorsAll{2};
end 

phaseNames = {'ITI', 'Cue', 'D1', 'Internal', 'D2', 'Speech'};
numPhases = length(phaseNames);

figBar = figure('units', 'normalized', 'outerposition', [0 0 0.3 0.4]);

for nMethod = 1:numMethods
    subplot(1, 2, nMethod)

    % Load data
    Data = load(fullfile(dataFolderPath, [unitRegion '_' regSub{nMethod} '_Speech_tuning.mat']));      
    randDotLocation = randi([98, 102], 1, Data.numSessions) * 0.01;
    
    % Determine sessions to include
    if strcmp(unitRegion, 'S1X') && strcmp(subjectId, 's3')
        sessionsToInclude = setdiff(1:Data.numSessions, 2); %session 2 S1 data not recorded
    else
        sessionsToInclude = 1:Data.numSessions;
    end
    
    numUnitsPerSession = Data.numUnitsPerSession;
    
    dataToPlot = cell2mat(Data.tuned_channels_per_phase(cueIdx, sessionsToInclude)') ./ numUnitsPerSession(cueIdx, sessionsToInclude)' * 100;
    
    sigVal = zeros(3, 3); 
    pValsAll = [0.05, 0.01, 0.001];

    pAll = [];
    % Perform t-tests between ITI-Cue, D1-Internal, D2-Speech
    n = 1;
    for phaseIdx = [1, 3, 5]
        [~, p_] = ttest(dataToPlot(:, phaseIdx), dataToPlot(:, phaseIdx + 1));
        sigVal(n, :) = p_ < pValsAll;
        pAll(n) = p_;
        cohens_d =utile.cohens_d(dataToPlot(:, phaseIdx), dataToPlot(:, phaseIdx + 1), 'paired');

        disp([regSub{nMethod} ':' phaseNames{phaseIdx} ' vs. ' phaseNames{phaseIdx+1}])
        disp(['p = ' num2str(p_,3) ' : cohens_d ' num2str(cohens_d,3)]);
      
        n = n+1;
    end
    [sigValITI, ITIpval] = ttest(dataToPlot(:, 1), dataToPlot(:, 4));
    % Display significant differences between ITI phase and Internal speech
    % phase
    if nnz(sigValITI)
        disp(['Using ' regSub{nMethod}])

        disp(['Internal phase is significantly different from ITI with p=' num2str(ITIpval) ' and cohens_d of ' num2str(utile.cohens_d(dataToPlot(:, 1), dataToPlot(:, 4), 'paired'),3)])
    end

    sigVal(isnan(sigVal)) = 0;
    sigValITI(isnan(sigValITI)) = 0;

    % Calculate error and confidence intervals
    errCi = utile.calculateErrCi(dataToPlot');

    % Plot bar graph with error bars
    dd = bar(mean(dataToPlot));
    hold on
    plot(1:numPhases, dataToPlot, 'Marker', 'o', 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'Color', [0, 0, 0], 'LineStyle', 'none');

    hold on 

    errorbar(1:numPhases, mean(dataToPlot), errCi(1, :), 'LineStyle', 'none', 'Color', [0 0 0])
    dd(1).FaceColor = plottingColors; 
    ylim([0 65])
    xticklabels(phaseNames)
    xtickangle(-45)
    ylabel('Tuned units [%]')
    title('Audio Cue')
    maxValA = mean(dataToPlot) + errCi(1, :);
    
    % Compute and display significant values and barlines  
    computeSigVal([1, 3, 5], sigVal, maxValA([2, 4, 6])); 
    set(gca, 'FontSize', 13)
    text(3, 60, 0, ['N = ' num2str(size(dataToPlot, 1))], 'FontSize', 15)
    title(regSub{nMethod})

end 

sgtitle([taskCue{1} ' cue - subject ' subjectId], 'fontSize', 18)

% Compute significant values and barlines  
function computeSigVal(phasesToCompute, sigValuesPerPhase, yData)

    for nPhaseToCompute = 1:length(phasesToCompute)
       nTmp = phasesToCompute(nPhaseToCompute);
       sigValuesTmp = sigValuesPerPhase(nPhaseToCompute, :);
       pSigns = {'*', '**', '***'};
       
       if nnz(sigValuesTmp) > 0
           pSign = pSigns(nnz(sigValuesTmp));
           text(nTmp + 0.5, yData(nPhaseToCompute) + 7, pSign, 'Color', 'black', 'LineWidth', 2, 'HorizontalAlignment', 'center');
           utile.sigline([nTmp, (nTmp + 1)], '', [], yData(nPhaseToCompute));
       end 
    end 
end
