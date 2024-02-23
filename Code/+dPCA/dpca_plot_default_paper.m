function dpca_plot_default_paper(data, time, yspan, explVar, compNum, events, signif, marg)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data, 'legend')
    
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return

    % if there is one parameter
    elseif length(time) == 3
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        hold on
        wordsToTest =  {'Battlefield','Cowboy','Python','Spoon','Swimming','Telephone','Bindip','Nifzig'};

        colorsCondition = utile.get_color_rgb_codes(wordsToTest);

        for l = 1:length(wordsToTest)
            hold on
            plot([0.5 1], [-(l+1) -(l+1)], 'color', colorsCondition{l}, 'LineWidth', 2)
            text(1.2, -(l+1), wordsToTest{l})
        end 
        
        %axis([0 3 -1 1.5+numOfStimuli])
        
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return

    % two parameters: stimulus and decision (decision can only have two
    % values)
    elseif length(time) == 4 && time(3) == 2
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end
        plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
        plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
        text(1.2, -2, 'Decision 1')
        text(1.2, -3, 'Decision 2')
        
        axis([0 3 -4.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
        
     elseif length(time) == 4 && time(3) == 5
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        StimuliNames = {'Image', 'Auditory', 'Written'};
        colors = cell2mat(utile.get_color_rgb_codes(StimuliNames)');
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, StimuliNames{f})
        end
        
        
         wordsToTest = {'Lateral', 'WritingTripod', 'MediumWrap',...
             'PalmarPinch','Sphere3Finger'};
        colorsCondition = utile.get_color_rgb_codes(wordsToTest);
        
        
        for l = 1:length(wordsToTest)
            hold on
            plot([0.5 1], [-(l+1) -(l+1)], 'color', colorsCondition{l}, 'LineWidth', 2)
            text(1.2, -(l+1), wordsToTest{l})
        end 
     
        
        %axis([0 3 -4.5 1.5+numOfStimuli])
        axis([0 3 -(1.5 + length(wordsToTest)) 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
        
        elseif   length(time) == 4 && time(3) == 8
            numOfStimuli = time(2); % time is used to pass size(data) for legend
            StimuliNames = {'Auditory', 'Written'};
            colors = cell2mat(utile.get_color_rgb_codes(StimuliNames)');
            hold on

            for f = 1:numOfStimuli
                plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
                text(1.2, f, StimuliNames{f})
            end


             wordsToTest =  {'Battlefield','Cowboy','Python','Spoon','Swimming','Telephone','Bindip','Nifzig'};

            colorsCondition = utile.get_color_rgb_codes(wordsToTest);

            for l = 1:length(wordsToTest)
                hold on
                plot([0.5 1], [-(l+1) -(l+1)], 'color', colorsCondition{l}, 'LineWidth', 2)
                text(1.2, -(l+1), wordsToTest{l})
            end 


            axis([0 3 -(1.5 + length(wordsToTest)) 1.5+numOfStimuli])
            set(gca, 'XTick', [])
            set(gca, 'YTick', [])
            set(gca,'Visible','off')
        return 
        
    % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(time)
    time = 1:size(data, ndims(data));
end
axis([time(1) time(end) yspan])
hold on

if ~isempty(explVar)
    title(['Component #' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'])
else
    title(['Component #' num2str(compNum)])
end

if  ~isempty(events) && length(events) == 6
    PhaseNames = {'ITI', 'Cue', 'D1', 'Imagined', 'D2', 'Speech'};
    for n_phase = 1:length(PhaseNames)
        xline(events(n_phase), 'k--',PhaseNames{n_phase}, 'LineWidth', 2, 'FontSize', 12)
    end 
    
elseif  ~isempty(events) && length(events) == 4
    PhaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
    for n_phase = 1:length(PhaseNames)
        xline(events(n_phase), 'k--',PhaseNames{n_phase}, 'LineWidth', 2, 'FontSize', 12)
    end 
elseif ~isempty(events)
    plot([events; events], yspan, 'Color', [0.6 0.6 0.6])

end


if ~isempty(signif)
    signif(signif==0) = nan;
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data) == 2
    % only time - plot it
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    wordsToTest =  {'Battlefield','Cowboy','Python','Spoon','Swimming','Telephone','Bindip','Nifzig'};
    colors = cell2mat(utile.get_color_rgb_codes(wordsToTest)');
            
    for f = 1:numOfStimuli
        plot(time, squeeze(data(1,f,:)), 'LineWidth', 2, 'Color', colors(f,:))   
    end 

elseif ndims(data) == 4 && size(data,3)==2
    % different stimuli in different colours and binary condition as
    % solid/dashed
    numOfStimuli = size(data, 2);
    colors = lines(numOfStimuli);

    for f=1:numOfStimuli 
        plot(time, squeeze(data(1, f, 1, :)), 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, :)), '--', 'color', colors(f,:), 'LineWidth', 2)
    end
    
elseif ndims(data) == 4 && size(data,3) == 5 && size(data,2) == 3

    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    
    
    dataImage = squeeze(data(:,1,:));
    dataAudio = squeeze(data(:,2,:));
    dataWritten = squeeze(data(:,3,:));

    if marg == 2 % = if its cue modality we want to plot different colors
         wordsToTest = {'Lateral', 'WritingTripod', 'MediumWrap',...
             'PalmarPinch','Sphere3Finger'};
        colorsCondition = utile.get_color_rgb_codes(wordsToTest);

        for i = 1:size(dataImage,2)
            plot(time, dataImage(:,i),'Color', colorsCondition{i},'LineWidth', 1.5);   
            hold on
            plot(time, dataAudio(:,i), 'Color', colorsCondition{i},'LineStyle', '--' ,'LineWidth', 1.5);   
            hold on
            plot(time, dataWritten(:,i), 'Color', colorsCondition{i},'LineStyle', '-.','LineWidth', 1.5);   
        end 
    else 

        
        wordsToTest = {'Lateral Image', 'Lateral Auditory' 'Lateral Written' ...
                 'WritingTripod Image', 'WritingTripod Auditory' 'WritingTripod Written' ...
                 'MediumWrap Image', 'MediumWrap Auditory' 'MediumWrap Written' ...
                 'PalmarPinch Image', 'PalmarPinch Auditory' 'PalmarPinch Written' ...
                 'Sphere3Finger Image', 'Sphere3Finger Auditory' 'Sphere3Finger Written'}; 
        colorsCondition = utile.get_color_rgb_codes(wordsToTest);
        for i = 1:size(dataImage,2)
            plot(time, dataImage(:,i),'Color', colorsCondition{3*i-2},'LineWidth', 1.5);   
            hold on
            plot(time, dataAudio(:,i), 'Color', colorsCondition{3*i-1},'LineWidth', 1.5);   
            hold on
            plot(time, dataWritten(:,i), 'Color', colorsCondition{:,3*i},'LineWidth', 1.5);   
        end 
    
    end 
    
    
elseif ndims(data) == 4 && size(data,3) == 4 && size(data,2) == 2

    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);

    dataGo = squeeze(data(:,1,:));
    dataNoGo = squeeze(data(:,2,:));

    if marg == 1
        wordsToTest = {'Go Trial', 'No Go Trial'};
        colorsCondition = utile.get_color_rgb_codes(wordsToTest);
        GrayColor = utile.get_color_rgb_codes({'redballoon'});
        for i = 1:size(dataGo,2)
            plot(time, dataGo(:,i),'Color', colorsCondition{1}, 'LineWidth', 1);
            hold on
            plot(time, dataNoGo(:,i), 'Color', colorsCondition{2},'LineStyle', '--', 'LineWidth', 1);   

            if i == 4
                plot(time, dataGo(:,i),'Color', GrayColor{1}, 'LineWidth', 1);
                hold on
                plot(time, dataNoGo(:,i), 'Color', GrayColor{1},'LineStyle', '--', 'LineWidth', 1);     
            end 
         end 
    else 

    wordsToTest = {'Lateral', 'WritingTripod', 'MediumWrap',...
         'redballoon'};
    colorsCondition = utile.get_color_rgb_codes(wordsToTest);

        for i = 1:size(dataGo,2)
            plot(time, dataGo(:,i),'Color', colorsCondition{i}, 'LineWidth', 1);
            hold on
            plot(time, dataNoGo(:,i), 'Color', colorsCondition{i},'LineStyle', '--' , 'LineWidth', 1);  
        end 
    end 
    
 elseif ndims(data) == 4 && size(data,2) == 2 && size(data,3) == 8
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    
    dataSoundCue = squeeze(data(:,1,:));
    dataWrittenCue = squeeze(data(:,2,:));
    
    wordsToTest =  {'Battlefield','Cowboy','Python','Spoon','Swimming','Telephone','Bindip','Nifzig'};

    colorsCondition = utile.get_color_rgb_codes(wordsToTest);

    if marg == 2
        for i = 1:size(dataSoundCue,2)
            plot(time, dataSoundCue(:,i),'Color', colorsCondition{i}, 'LineWidth', 1);
            hold on
            plot(time, dataWrittenCue(:,i), 'Color', colorsCondition{i},'LineStyle', '--' , 'LineWidth', 1);  
        end 
     
    else
        
        colorsCondition =utile.get_color_rgb_codes({'Auditory', 'Written'});
        
        for i = 1:size(dataSoundCue,2)
            plot(time, dataSoundCue(:,i),'Color', colorsCondition{1}, 'LineWidth', 1);
            hold on
            plot(time, dataWrittenCue(:,i), 'Color', colorsCondition{2},'LineStyle', '--' , 'LineWidth', 1);  
        end 
        
    end 
    

else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    
    plot(time, data, 'LineWidth', 2)    
end
