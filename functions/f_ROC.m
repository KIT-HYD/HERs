function [TP, TN, FP, FN, TPR, TNR, FPR, Accuracy, F1_score, MCC, DistanceROC] = f_ROC(true_z, predicted)
%ROC curve

    TP = sum(true_z .* predicted);        %true positive
    TN = sum((1-true_z) .* (1-predicted));  %true negatve
    FP = sum((1-true_z) .* predicted);      %false positive
    FN = sum((true_z) .* (1-predicted));    %false negatve
    
    P = (TP + FN);                      %positive
    N = (TN + FP);                      %negative
    
    TPR = TP/P;                   %true positive rate or sensitivity
    TNR = TN/N;                   %true negative rate or specificity
    FPR = FP/(FP+TN); %1-TNR      %false positive rate
  
    Accuracy = (TP + TN) / (P + N);
    F1_score = 2*TP /(2*TP + FP + FN);
    MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); %Matthews correlation coefficient
    DistanceROC = sqrt((1-TPR)^2 + (0-FPR)^2);
    
end
