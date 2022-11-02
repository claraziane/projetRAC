%% This function extracts participants' cadence from step onset time values
%
% Input variables:
% -stepOnsetAll: vector (length = number of steps) of left and right step onset time values (in frames)
% -Freq: acquisition frequency
%
% Output variables:
% -cadence: vector (length = number of right steps) of right step onset time values (in frames)
% -stepFreq: single value for step frequency
%
% C. Ziane

function[cadence, stepFreq] = getCadence(stepOnsetAll, Freq)

nSteps   = length(stepOnsetAll);
stepFreq = nSteps / ((stepOnsetAll(end) - stepOnsetAll(1))/Freq);
cadence  = stepFreq * 60;

end