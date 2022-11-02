%% Randomize condition order
clear;
close;
clc;

Participant = input('What is the participant''s code name ?', 's');
if ~exist(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant], 'dir')
    mkdir(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant])
end

conditionsOrder = {'noRAC_preferredRest'; []; []; []; 'RAC_preferredRest';  [];   [];   []};
Conditions = {'preferredWalk'; 'slowWalk'; 'fastWalk'};

conditionsRandom = randperm(3);
i = 1;

for iCond = 1:length(conditionsOrder)

    if iCond > 1 && iCond < 5
        conditionsOrder(iCond,1) = strcat('noRAC_', Conditions(conditionsRandom(i)));
        i = i+1;
    elseif iCond > 5
        conditionsOrder(iCond,1) = strcat('RAC_', Conditions(conditionsRandom(i)));
        i = i+1;
    end

    if iCond == 4
        conditionsRandom = randperm(3);
        i = 1;
    end
end

save(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant '/conditionsOrder.mat'], 'conditionsOrder');