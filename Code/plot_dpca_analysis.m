

%% compute and plot cross-classification decoding results (Figure 6A,B of manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); %add folder to search path 
clc
clear all
%close all
 
%% parameters
subject_id= 's3';
unit_region = 'SMG'; %SMG or S1

dataFolderPath = fullfile(pwd,'Database', subject_id,'Data');

tableName = ['Table_' subject_id  '_Speech.mat'];
Data = load(fullfile(dataFolderPath,tableName)); 
Go_data = Data.Go_data; 

brainRegionIdx = ismember(Go_data.dPCA{6}, unit_region);
firingRates = Go_data.dPCA{brainRegionIdx};

N = size(firingRates,1);   % number of neurons
S = size(firingRates,2);   % number of cues: auditory cue and written cue
D = size(firingRates,3);   % number of spoken words
T = size(firingRates,4);   % number of time points
E = size(firingRates,5);   % maximal number of trial repetitions

trialNum = E*ones(N,S,D); %number of trials for each neuron in each S,D condition

firingRatesAverage = nanmean(firingRates,5);

time = (1:T) / 20; % neural data recorded in 50ms time bins

ifSimultaneousRecording = true;  % change this to simulate simultaneous 
                                 % recordings (they imply the same number 
                                 % of trials for each neuron)

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines

timePhaseLabels = Go_data.time_phase_labels{1};
timeIdx(1) = 1; 
timeIdx(1,2:length(unique(timePhaseLabels))) = (find(diff(timePhaseLabels))+1)';
timeEvents = time(timeIdx);

%% Define parameter grouping

%Replace nan values by 0 (dPCA not compatible with nan values). 
warning('removing nan values')
firingRates(isnan(firingRates)) = 0;

% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
        end
    end
end


% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

%combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};

if strcmp(subject_id, 's2')
    margNames = {'Cue Modality', 'Word', 'Timing', 'Cue/Word Interaction'};
    combinedParams = {{1, [1 3]}, {2, [2 3]},{3}, {[1 2],[1 2 3]}};

elseif  strcmp(subject_id, 's3')
    margNames = {'Word', 'Word/Timing', 'Timing'};
    combinedParams = {{1}, {[1 2]}, {2}};  
    firingRatesAverage = squeeze(firingRatesAverage);
    firingRates = squeeze(firingRates);
    trialNum = squeeze(trialNum);

end 

margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
%timeEvents = time(round(length(time)/2));



%% Step 4: dPCA with regularization

%compute or load optimal lambda values 
lambdaFile = fullfile(pwd , 'Code', '+dPCA','optimalLambdas',  [subject_id '_' unit_region '_optimalLambdas.mat']);
if ~exist(lambdaFile)
    
    optimalLambda = dPCA.dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 10, ...  
        'filename', ['optimalLambdas/' subject_id '_' unit_region '_optimalLambdas.mat']);

else
    optimalLambda = load(lambdaFile, 'optimalLambda');
end 

Cnoise = dPCA.dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dPCA.dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dPCA.dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dPCA.dpca_plot(firingRatesAverage, W, V, @dPCA.dpca_plot_default_paper, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

