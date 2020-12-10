function [accuracy] = CalcAccuracy(Z_model,Z_actual)
%CALCACCURACY Calculate Accuracy
%   Standard accuracy calculation
diff = abs(Z_model - Z_actual);
ratio = diff./Z_actual;
accuracy = max(1 - ratio,0);
accuracy(isnan(accuracy)) = 0;
accuracy = mean(accuracy,'all');
end

