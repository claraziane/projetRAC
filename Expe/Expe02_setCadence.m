%% Determine cadence for all conditions
% clear;
% close;
% clc;

Participant = input('What is the participant''s code name ?', 's');
if ~exist(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant], 'dir')
    mkdir(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant])
end

Cadence.preferredWalk = input('What is the participant''s preferred cadence (steps/min) ?');
Speed.  preferredWalk = input('What is the treadmill speed at preferred cadence (m/s) ?');

Cadence.slowWalk = round(Cadence.preferredWalk * 0.9);
Cadence.fastWalk = round(Cadence.preferredWalk * 1.1);

Speed.fastWalk = input('What is the treadmill speed at fast cadence (m/s) ?');
Speed.slowWalk = input('What is the treadmill speed at slow cadence (m/s) ?');

save(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/cprojetRAC/DATA/Behavioural/' Participant '/Cadence&Speed.mat'], 'Speed', 'Cadence');