%calculate standard error of the mean

function err_ci = calculateErrCi(Data)
    NA = size(Data,2);
    SEM_A = std(Data') / sqrt(NA); % Standard Error Of The Mean
    CI95 = tinv([0.025 0.975], NA-1);   
    err_ci =  abs(bsxfun(@times, SEM_A, CI95(:)));
end 
