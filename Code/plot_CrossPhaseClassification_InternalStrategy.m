%% compute and plot cross-classification decoding results (Figure 6A,B of manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); %add folder to search path 
clc
clear all
%close all

%% parameters 
subjectId= 's2'; % only for s2
unitRegion = 'SMG'; % only for SMG

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd, 'Database', subjectId, 'Data', 'ClassificationData', 'CrossPhaseInternalStrategy');     
cueforSpeechAll = {'Audio_Audio','Audio_Written', 'Written_Written','Written_Audio'};
InternalStrategy = {'Sound','Visual', 'Visual', 'Sound'};
subset_idx = [1,2,4,6];
errAll = zeros(4,4); 
for n_strategies = 1:length(cueforSpeechAll)
	cueforSpeech = cueforSpeechAll{n_strategies};
    errTmp = load(fullfile(dataFolderPath, cueforSpeech));
    shuffleTmp = load(fullfile(dataFolderPath,[cueforSpeech 'Shuffle']));

    errAll(n_strategies,:) = errTmp.errTesting(4,subset_idx);
    shuffleAll(:,n_strategies,:) = shuffleTmp.errTesting(:,4,subset_idx); 
end

prc_vals = [97,99.5,99.95];
sigvals = {};
for n_prc = 1:length(prc_vals)
    sigvals{n_prc} = squeeze(prctile(shuffleAll,prc_vals(n_prc))) < errAll; 
end 


%% Only plot ITI Cue InternalSpeech and Speech
PhaseNames = {'ITI', 'Cue', 'Internal', 'Speech'};
ColorForPlot =utile.get_color_rgb_codes(PhaseNames);

fig2 = figure('units','normalized','outerposition',[0 0 0.3, 0.4]);

subplot(1,2,1)

ba = bar(errAll(1:2,:), 'FaceColor', 'flat');
for ck = 1:length(ColorForPlot); ba(ck).CData = ColorForPlot{ck}; end
ylabel('Classification accuracy')
xticklabels(InternalStrategy(1:2));
xlabel('Internal Strategy')
title('Audio cue')
set(gca, 'FontSize', 12)

subplot(1,2,2)

ba = bar(errAll(3:4,:), 'FaceColor', 'flat');
for ck = 1:length(ColorForPlot); ba(ck).CData = ColorForPlot{ck}; end
xticklabels(InternalStrategy(3:4));
xlabel('Internal Strategy')
title('Written cue')
set(gca, 'ytick', [])
set(gca, 'FontSize', 12)

    




