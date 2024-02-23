function [cohen_d] = cohens_d(condition1, condition2, testType)
    
    %return absolute value of paired conhen's d test for either paired or
    %independent sets
        
    mean_condition1 = nanmean(condition1);
    mean_condition2 = nanmean(condition2);
    std_condition1 = nanstd(condition1);
    std_condition2 = nanstd(condition2);
    
    n1 = length(condition1);
    n2 = length(condition2);
    
    
    switch testType
    
        case 'independent'
           
            % Compute pooled standard deviation
            pooled_std = sqrt(((n1 - 1) * std_condition1^2 + (n2 - 1) * std_condition2^2) / (n1 + n2 - 2));
        
            % Compute Cohen's d
            cohen_d = abs(mean_condition1 - mean_condition2) / pooled_std;
    
        case 'paired'
            % compute differences of conditions
            condition_delta   = condition1 - condition2;     
            % compute standard deviation of deltas
            std_delta = nanstd(condition_delta);              
            mean_diff = mean(condition_delta);
            cohen_d        =  abs(mean_diff / std_delta);   % Cohen's d (paired version)
    
        otherwise
           error('Provide testType - independent or paired')
    
    
    end
      
end
