
%% compute and plot cross-classification decoding results (Figure 6A,B of manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

addpath(genpath(pwd)); %add folder to search path 
clc
clear all
%close all

%% parameters 
subjectId= 's2'; % only for s2
unitRegion = 'SMG'; % only for SMG
%taskCue = 'Auditory'; % Written or Auditory 
taskCue = 'Written';

% Define the data folder path based on subject ID
dataFolderPath = fullfile(pwd,'Database', subjectId,'Data','ClassificationData', 'CrossPhaseClassification');

%% plot data 



ClassificationData = load(fullfile(dataFolderPath,[unitRegion '_' taskCue '_Speech_Data']));

TestingErrorOrig = ClassificationData.errTest;
% plot all phases
numPhases = size(TestingErrorOrig,4);

%TestingErrorOrig
keepMeanAcc = squeeze(TestingErrorOrig);
dataPhasesToPlot = [1,2,4,6];

PhaseNames = {'ITI', 'Cue','D1', 'Internal','D2', 'Speech'};
PhaseNamesSubset = PhaseNames(dataPhasesToPlot);

ColorModel = utile.get_color_rgb_codes(PhaseNamesSubset);
numPhasesToPlot = numel(dataPhasesToPlot); 

keepMeanAcc = (1- keepMeanAcc([dataPhasesToPlot],[dataPhasesToPlot],:))*100;

modelError = [];
maxVal = [];

for nPhase = 1:numPhasesToPlot
    err_ci =utile.calculateErrCi(squeeze(keepMeanAcc(nPhase,:,:)));
    modelError(nPhase,:) = err_ci(1,:);
    maxVal(nPhase,:) = mean(squeeze(keepMeanAcc(nPhase,:,:)),2) + err_ci(1,:)' ;

end 

sigVal = cell(1,numPhasesToPlot); 
pVal = cell(1,numPhasesToPlot); 
cohensD = cell(1,numPhasesToPlot); 
%compute ttests between phases
for nSubPhase = 1:numel(dataPhasesToPlot)
    
    sigVal{nSubPhase} = zeros(numPhasesToPlot,numPhasesToPlot); 
    pVal{nSubPhase} = ones(numPhasesToPlot,numPhasesToPlot); 
    cohensD{nSubPhase} = ones(numPhasesToPlot,numPhasesToPlot);  

    for nPhase = 1:numPhasesToPlot
        
        for tPhase = setdiff(1:numPhasesToPlot,nPhase)
            meanTmp = squeeze(keepMeanAcc(nSubPhase,:,:));           
            [~,pVal{nSubPhase}(nPhase,tPhase)] = ttest(meanTmp(nPhase,:), meanTmp(tPhase,:));
            cohensD{nSubPhase}(nPhase,tPhase) = utile.cohens_d(meanTmp(nPhase,:), meanTmp(tPhase,:), 'paired');
         end 
    end 
    
    maskU = triu(true(size(pVal{nSubPhase})),+1);

    % correct p-values for multiple comparisons
    p_vals = pVal{nSubPhase}(maskU);
    [a,~,p_adj] = utile.MultipleComparisonsCorrection(p_vals,'method', 'fdr'); 
    
    pVal{nSubPhase}(maskU) = p_adj;
    pVal{nSubPhase}(~maskU) = inf;
    cohensD{nSubPhase}(~maskU) = inf;

    maskL = tril(true(size(pVal{nSubPhase})),-1);
     
    sigVal{nSubPhase} = pVal{nSubPhase} < 0.05;
    % lower bound = repetition, put to zero
    sigVal{nSubPhase}(maskL) = 0;
    % first row, ITI, put to zero to not plot significance between ITI and other task phases
    sigVal{nSubPhase}(1,:) = 0;  

end 

modelSeries = mean(keepMeanAcc,3); 
fig = figure();
b = bar(modelSeries, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(modelSeries);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for nPhase = 1:nbars
    x(nPhase,:) = b(nPhase).XEndPoints;
    b(nPhase).FaceColor = ColorModel{nPhase};
end
% Plot the errorbars
x = x';
z = errorbar(x,modelSeries,modelError,'k','linestyle','none');
hold off
title([' Cross Phase Classification - ' taskCue])   
xticks(1:numPhases);
xticklabels(PhaseNamesSubset)
xtickangle(-45)
ylabel('Classification accuracy [%]');
xlabel('Training phase');

set(gca, 'FontSize', 15)
ylim([0 100])
add = 0;
% plot bar figures with sig. lines

for nPhase = 1:numPhasesToPlot

    sigValTmp = sigVal{nPhase};
    pValTmp = pVal{nPhase};
    cohensDTmp = cohensD{nPhase};

    for mPhase = 1:numPhasesToPlot
        
        toPlot = find(sigValTmp(mPhase,:)); 
        toPlotP = pValTmp(mPhase,:); 
        cohensDToPlot = cohensDTmp(mPhase,:);
        for tPhase = toPlot  
            maxV = max(maxVal(nPhase,tPhase), maxVal(nPhase,mPhase));
    
            utile.sigline([x(nPhase,mPhase),x(nPhase,tPhase)],'', [], maxV + add);
            text((x(nPhase,mPhase)+ x(nPhase,tPhase)) /2,maxV +add +2,plot_p_val(toPlotP(tPhase)),'FontSize',15, 'HorizontalAlignment', 'center');
            disp([ ' Training phase ' PhaseNamesSubset{nPhase} ' : ' PhaseNamesSubset{mPhase}  ' vs. ' PhaseNamesSubset{tPhase} ' p = ' num2str(toPlotP(tPhase),3) ' - cohens_d ' num2str(cohensDToPlot(tPhase),3)])

            add = add+7;
        end 
        
    end 

    add = 0;
    
end 

%% Helper function

function p_val_sign = plot_p_val(p_val)  
    if p_val < 0.001
        p_val_sign = '***';
    elseif p_val < 0.01
        p_val_sign = '**';          
    elseif p_val < 0.05
        p_val_sign = '*';
    end
end 
