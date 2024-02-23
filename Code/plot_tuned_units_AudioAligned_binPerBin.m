%% plot regression tuning in 50ms time bins (Figure 3A,D, S3A,B in manuscript)

%% Important: run code while being in folder 'DatabaseInternalSpeech' or change dataFolder pathway

% Add current folder and its subfolders to the search path
addpath(genpath(pwd)); 
clc;
clear all;

subject_id = 's2'; % change to participant s3
unit_region = 'SMG'; 
flagRegressionTuning = true; %True = Regression. False = Kruskwal wallis

dataFolderPath = fullfile(pwd, 'Database', subject_id, 'Data', 'AlignedTuningData');

%chose cue type:
if strcmp(subject_id, 's2')
   taskCuesAll = {'Auditory', 'Written'};
   plottingColors = utile.get_color_rgb_codes({'Speech', 'Written'});
   cueOrder = [2,1];
elseif strcmp(subject_id, 's3')
    taskCuesAll = {'Written'}; %only Written cue trials for s3
    plottingColors = utile.get_color_rgb_codes({'Written'});
    cueOrder = [1];
else
    keyboard % unknown subject
end

numCues = numel(taskCuesAll);

reg_sub = {'KruskalWallis','Regression'};
%% plot all sessions with 95% confidence interval
    
fontSize =15; 

Data = load(fullfile(dataFolderPath,[unit_region '_' reg_sub{flagRegressionTuning + 1}  '_Speech_tuning.mat' ]));        

numUnitsPerSession = Data.numUnitsPerSession;
h = [];
Data.numSessions = size(Data.numUnitsPerSession,2);
if strcmp(unit_region, 'S1X') && strcmp(subject_id, 's3')
    sessionsToInclude = setdiff(1:Data.numSessions,2);
else
    sessionsToInclude = 1:Data.numSessions;
end 

binPerbinPerPhase= figure('units','normalized','outerposition',[0 0 0.35 0.43]);
timePhaseLabels = median(reshape(Data.timePhaseLabelsPerSession, Data.numCues*Data.numSessions, []));
phaseNames ={'ITI', 'Cue','D1','Internal','D2', 'Speech'};
numPhases = length(phaseNames);

for n_phase = 1:numPhases 
    
    if n_phase  == 1
        subplot(1,15,1:4)
    elseif n_phase == 2
        subplot(1,15,5:7)
    elseif n_phase == 3
        subplot(1,15,8)
    elseif n_phase == 4
        subplot(1,15,9:11)
    elseif n_phase == 5
        subplot(1,15,12)
    elseif n_phase == 6
        subplot(1,15,13:15)
    end 
     
  
    for i = 1:numCues
        special_color = plottingColors{i};

        DataToPlot = (cell2mat(Data.sum_bin_all(i,sessionsToInclude))./numUnitsPerSession(i,sessionsToInclude))';

        if i == 1 && n_phase == 2
            %relace zeros with nans for plotting auditory trials
            DataToPlot(:,61:73) = nan;
        end 
        yLim = 60; 
      
        DataToPlot = DataToPlot(:,timePhaseLabels == n_phase);
        %95 Percent Confidence interval 
        x = 1:length(DataToPlot);                                  
        acc = DataToPlot*100;
        NA = size(acc,1);
        meanAcc = mean(acc);
        SEM_A = std(acc) / sqrt(NA); % Standard Error Of The Mean
        CI95_A = SEM_A * tinv(0.975, NA-1);  
        time_idx = (1:length(timePhaseLabels))*0.05;
        time_idx = time_idx(timePhaseLabels == n_phase);
    
        ER =utile.shadedErrorBar(time_idx,meanAcc, CI95_A);
        ER.mainLine.Color = special_color;
        ER.patch.FaceColor = special_color;
        ER.edge(1).Color = special_color;
        ER.edge(2).Color = special_color;  
        hold on
        h{i} = plot(time_idx,meanAcc, 'Linewidth', 3);
        h{i}.Color = special_color;
        hold on       
        ylim([0 yLim])
        xlim([time_idx(1), time_idx(end) ])
        
    end 
    
    uistack(h{1},'top')
    
    if n_phase == 1
        ylabel('Tuned units [%]', 'FontSize', 13)
    end 
     
    if n_phase ~=1
        set(gca,'YTick', [])
    end 
    
    ax = gca;
    ax.XTick = unique( round(ax.XTick) );
    title(phaseNames{n_phase});        
    set(gca, 'Fontsize', fontSize)

end 

taskCuesAllLabel = cellfun(@(x) [x ' - 95% CI.'], taskCuesAll, 'UniformOutput', false);
legend([h{:}], taskCuesAllLabel, 'FontSize', fontSize)
CondToTest = utile.image2class_simple(unique(Data.labelsPerSession));
    
