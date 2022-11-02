%% Calculating the stability index
% 1. Compute RESS component :
%   a. sequence EEG in trials according to beat onsets (-100:500 ms)
%   b. compute covariance matrices S (from narrow-filtered data) and R (from broad-band signal)
%   c. perform Generalized Eigen Decomposition (GED)
% 2. Compute stability index
%   a. Filter RESS component
%   b. Transform filtered RESS into analytical signal using Hilbert transform
%   c. Compute phase angles
%   d. Extract instantaneous frequencies
%   e. Compute stability index (standard deviation of instantaneous frequencies)

clear;
close all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
%     pathData    = '\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\';
%     pathResults = '\\10.89.24.15\j\Projet_RAC\dataAnalysis\Results\EEG\';
    addpath('C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0');
    addpath('C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\dataAnalysis\Toolbox\GED-master');
    addpath '\\10.89.24.15\j\Projet_RAC\dataAnalysis\Scripts\Clara\EEG\Functions'
    addpath '\\10.89.24.15\j\Projet_RAC\dataAnalysis\Scripts\Clara\Stats\Functions'
    load([pathData 'RAC\All\RAC.mat']);
    load([pathData 'Kinetics\All\stepData.mat']);
else
    pathData    = '/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/'; %j
    pathEEG     = '/Volumes/Projet_RAC/DATA/Preprocessed/EEG/'; %q
    pathResults = '/Volumes/Seagate/Clara/UdeM/projetRAC/dataAnalysis/Results/EEG/';
%     pathResults = '/Volumes/Clara/UdeM/S2M/projetRAC/dataAnalysis/Results/EEG/';
    addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1');
    addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/GED-master');
    addpath '/Volumes/10.89.24.15/Projet_RAC/dataAnalysis/Scripts/Clara/EEG/Functions'
    %% 
    load([pathData 'RAC/All/RAC.mat']);
    load([pathData 'Kinetics/All/stepData.mat']);
    load([pathResults 'All/flanking_noEpoch_icaClean/resultsEEG.mat']);
%     load('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/RESS-main/Data/lf.mat')
end

Participants = {'P02'; 'P03';'P04'; 'P05'; 'P10'; 'P11'};
Conditions = {'RAC_preferredWalk';   'RAC_slowWalk';  'RAC_fastWalk';    'RAC_preferredRest'; ...
            'noRAC_preferredWalk'; 'noRAC_slowWalk'; 'noRAC_fastWalk'; 'noRAC_preferredRest'};
extensionRoot = '_chanClean.set';
analysisType  = 'flanking_noEpoch/';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 1:length(Participants)

    for iCondition = 1:length(Conditions)

        %% Load EEG

        directory = fullfile(pathEEG, Participants{iParticipant}, '/');
        directoryResults = fullfile(pathResults, Participants{iParticipant}, '/', Conditions{iCondition}, '/', analysisType);
        
        % Create folder in participant's result folder if does not exist
        if ~exist(directoryResults, 'dir')
            mkdir(directoryResults)
        end

        if ~exist(fullfile(directory, strcat(Conditions{iCondition}, extensionRoot)), 'file')
            directory = fullfile(pathData, 'EEG/', Participants{iParticipant}, '/');
        end
        EEG = pop_loadset('filename', strcat(Conditions{iCondition}, extensionRoot), 'filepath', directory); % Loads an EEG datas
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','on'); % Edits/saves EEG dataset structure information

        % Resample
        EEG = pop_resample(EEG, EEG.srate/3);
        EEG = eeg_checkset(EEG);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        dataRaw  = double(EEG.data);

%         % Epoch with auditory cues
%         EEG = pop_epoch(EEG, {'RAC'}, [-0.5 1], 'epochinfo', 'yes');
%         EEG = eeg_checkset(EEG);
%         EEG = pop_rmbase(EEG, [-500 0] ,[]);
%         EEG = eeg_checkset(EEG);
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','on');

        % Increase data precision (from single to double) for numerical stability
        data = double(EEG.data);

        %% specify parameters

        if strcmp(Conditions{iCondition}(1:5), 'noRAC')
            if strcmp(Conditions{iCondition}(end-3:end), 'Walk')
                stepFreq = round(stepData.([Participants{iParticipant}]).([Conditions{iCondition}]).stepFreq,2);
                sFreq    = stepFreq;
            elseif strcmp(Conditions{iCondition}(end-3:end), 'Rest')
                sFreq = round(RAC.([Participants{iParticipant}]).RAC_preferredRest.beatFrequency,2);
            end
        elseif strcmp(Conditions{iCondition}(1:3), 'RAC')
            stimFreq = round(RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).beatFrequency,2);
            if strcmp(Conditions{iCondition}(end-3:end), 'Walk')
                stepFreq = round(stepData.([Participants{iParticipant}]).([Conditions{iCondition}]).stepFreq,2);
                sFreq    = stepFreq;
            elseif strcmp(Conditions{iCondition}(end-3:end), 'Rest')
                sFreq = stimFreq;
            end
        end

        % Electrode used for 'best-electrode' analyses
        electrode = 'cz';

        % RESS parameters
        loFreq = 1; % distance of neighboring frequencies away from peak frequency, +/- in Hz
        hiFreq = 1;
        sFWHM  = 0.5; % FWHM at peak frequency
        loFWHM = 0.5; % FWHM of the neighboring frequencies
        hiFWHM = 0.5;

        % FFT parameters
        nFFT    = ceil( EEG.srate/.01 ); % .1 Hz resolution
        timeWin = dsearchn(EEG.times',[10000 299000]'); % we use data from .5-10 seconds, to avoid using the stimulus transient

        % FFT
        dataX = mean(abs(fft(data(:,timeWin(1):timeWin(2),:),nFFT,2)/diff(timeWin) ).^2,3);
        Hz    = linspace(0,EEG.srate,nFFT);

        %% Compute covariance matrices

        % compute covariance matrix at peak frequency
        dataS = filterFGx(data,EEG.srate,sFreq,sFWHM);
%         dataS = reshape( dataS(:,tidx(1):tidx(2),:), EEG.nbchan,[] );
        dataS = bsxfun(@minus,dataS,mean(dataS,2));
        covarianceS  = (dataS*dataS')/diff(timeWin);

        % compute covariance matrix for lower neighbor
        dataLo = filterFGx(data,EEG.srate,sFreq-loFreq,loFWHM);
%         dataLo = reshape( dataLo(:,tidx(1):tidx(2),:), EEG.nbchan,[] );
        dataLo = bsxfun(@minus,dataLo,mean(dataLo,2));
        covarianceLo  = (dataLo*dataLo')/diff(timeWin);

        % compute covariance matrix for upper neighbor
        dataHi = filterFGx(data,EEG.srate,sFreq+hiFreq,hiFWHM);
%         dataHi = reshape( dataHi(:,tidx(1):tidx(2),:), EEG.nbchan,[] );
        dataHi = bsxfun(@minus,dataHi,mean(dataHi,2));
        covarianceHi  = (dataHi*dataHi')/diff(timeWin);

        covarianceR = (covarianceHi+covarianceLo)/2;

        % Apply regularization to R
        regulFactor = .01;
        covarianceR = (1-regulFactor)*covarianceR + regulFactor*mean(eig(covarianceR))*eye(size(covarianceR));

%         clim = [-1 1]*10;
%         figure(3), clf
%         subplot(221); imagesc(covarianceS);
%         axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%         title('Covariance Matrix S');
%         subplot(222); imagesc(covarianceR);
%         axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%         title('Covariance Matrix R');
%         subplot(223); imagesc(covarianceLo);
%         axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%         title('Covariance Matrix Low');
%         subplot(224); imagesc(covarianceHi);
%         axis square; set(gca,'clim',clim); xlabel('Channels'), ylabel('Channels'); colorbar
%         title('Covariance Matrix Hi');
%         saveas(figure(3), [directoryResults '/fig_covarianceMatrices.png']);

        % Generalized eigendecomposition
        [W,Lambdas] = eig(covarianceS, covarianceR);
        [lambdaSorted, lambdaIndex] = sort(diag(Lambdas), 'descend');
        [~,compMax] = max(diag(Lambdas)); % find maximum component
%         compMax = compMax-1;
        W = bsxfun(@rdivide, W, sqrt(sum(W.^2,1))); % normalize vectors (not really necessary, but OK)

        % Extract components and force sign
        maps = covarianceS * W / (W' * covarianceS * W);

        if strcmp(Conditions{iCondition}(end-3:end), 'Walk')
            % Force Cz to be positive
           elecSign = sign(maps(strcmpi(electrode,{EEG.chanlocs.labels}),compMax));
            maps = maps * elecSign;
%             elecSign = -1;
        elseif strcmp(Conditions{iCondition}(end-3:end), 'Rest')
            % Force electrode with highest activity to be positive
            [~,idx] = max(abs(maps(:,compMax))); % Find biggest component
            maps = maps * sign(maps(idx,compMax)); % Force to positive sign
        end

        %% Reconstruct RESS component time series

        RESS = zeros(EEG.pnts, size(data,3));
        for iTrial = 1:size(data,3)
            RESS(:,iTrial) = W(:,compMax)'*squeeze(data(:,:,iTrial));
        end
        RESStime = W(:,compMax)'*dataRaw;

        [f, FFT] = spec_fft(RESStime, EEG.srate, 0);
        [M, fIndex] = max(f);
        timeVector = linspace(1, length(RESStime)/EEG.srate, length(RESStime));

%         figure(4)
%         subplot(221); map2plot = maps(:,compMax); topoplot(map2plot./max(map2plot), EEG.chanlocs, 'maplimits', [-0.7 0.7], 'numcontour',0,'conv','off','electrodes','on','shading','interp'); colorbar;
%         title([ 'RESS for ' num2str(sFreq) ' Hz' ], 'FontSize', 14);
%         subplot(222); plot(timeVector, RESStime);
%         set(gca, 'xlim', [timeVector(1) timeVector(end)]);
%         xlabel({'Time (s)'}, 'FontSize', 14),
%         title([ 'RESS for ' num2str(sFreq) ' Hz in the time domain' ], 'FontSize', 14);
%         subplot(2,2,[3:4]); plot(FFT,f); xlim = [0 25]; set(gca,'xlim',xlim); xlabel('Frequency (Hz)', 'FontSize', 14) ; ylabel('Power', 'FontSize', 14);...
%             legend(['Peak frequency = ' num2str(FFT(fIndex))], 'FontSize', 14);
%         title('FFT of RESS Component', 'FontSize', 14);
%         saveas(figure(4), [directoryResults '/fig_ressTopo.png']);

        i = 1;
%         figure(5), clf
%         subplot(211); plot(lambdaSorted,'ks-','markersize',10,'markerfacecolor','w');
%         xlabel('Component'); ylabel('\lambda');
        RESSsorted = zeros(5,size(dataRaw,2));
        for iComp = [1:4 length(Lambdas)]
            RESSsorted(i,:) = W(:,lambdaIndex(iComp))'*dataRaw;

%             subplot(2,5,5+i); map2plot = maps(:,lambdaIndex(iComp)); topoplotIndie(map2plot./max(map2plot), EEG.chanlocs, 'maplimits', [-1 1], 'numcontour',0,'conv','on','electrodes','off','shading','interp');
%             title([ 'Component ' num2str(iComp) ])
            i = i+1;
        end
%         colormap jet
%         saveas(figure(5), [directoryResults '/fig_ressComponents.png']);

        %% Compute SNR spectrum
        RESSx    = mean(abs(fft(RESS(timeWin(1):timeWin(2),:),nFFT,1)/diff(timeWin)).^2,2);
        elecx    = dataX(strcmpi(electrode,{EEG.chanlocs.labels}),:,:);

        [snrR,snrE] = deal(zeros(size(Hz)));
        skipbins =  5;
        numbins  = 20+skipbins;

        % Loop over frequencies and compute SNR
        for iHz = numbins+1:length(Hz)-numbins-1
            numer = RESSx(iHz);
            denom = mean(RESSx([iHz-numbins:iHz-skipbins iHz+skipbins:iHz+numbins]));
            snrR(iHz) = numer./denom;

            numer = elecx(iHz);
            denom = mean( elecx([iHz-numbins:iHz-skipbins iHz+skipbins:iHz+numbins]) );
            snrE(iHz) = numer./denom;
        end

%         figure(6), clf
%         xlim = [0.5 10.5];
%         subplot(2,2,1); map2plot = maps(:,compMax); topoplotIndie(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
%         title([ 'RESS for ' num2str(sFreq) ' Hz' ]);
%         subplot(2,2,2); map2plot = dataX(:,dsearchn(Hz',sFreq)); topoplotIndie(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','emarker2',{find(strcmpi({EEG.chanlocs.labels},electrode)) 'o' 'w' 4},'shading','interp');
%         title([ 'RESS for ' num2str(sFreq) ' Hz' ]);
%         title([ 'Electrode power at ' num2str(sFreq) ' Hz' ]);
%         subplot(2,2,[3:4]); plot(Hz,snrR,'ro-','linew',1,'markersize',5,'markerface','w'); hold on;
%         plot(Hz,snrE,'ko-','linew',1,'markersize',5,'markerface','w');
%         set(gca,'xlim',xlim); xlabel('Frequency (Hz)'), ylabel('SNR'); legend({'RESS';electrode}); clear xlim
%         saveas(figure(6), [directoryResults '/fig_ressVSelectrode.png']);

        %% Computing the stability index

        % Gaussian filering RESS component
        ressFiltered = filterFGx(RESStime, EEG.srate, sFreq, sFWHM);

        % Compute Hilbert Transform
        ressHilbert = hilbert(ressFiltered);

        % Extract phase angles
        ressPhaseAngle = angle(ressHilbert);
        ressPhaseAngle = unwrap(ressPhaseAngle);
        ressPhaseHz    = (EEG.srate*diff(ressPhaseAngle)) / (2*pi);

        % Apply a sliding moving median with a window width of 400 samples
        nOrder = 10;
        orders = linspace(10,400,nOrder)/2;
        orders = round(orders/(1000/EEG.srate));
        phasedMed = zeros(length(orders), length(ressPhaseHz));
        for iOrder = 1:nOrder
            for iTime = 1:length(ressPhaseHz)
                temp = sort(ressPhaseHz(max(iTime-orders(iOrder),1):min(iTime+orders(iOrder),length(ressPhaseHz)-1)));
                phasedMed(iOrder,iTime) = temp(floor(numel(temp)/2)+1);
            end
        end
        phaseMedFilt = mean(phasedMed);

        % Compute stability index
        stabilityIndex = std(phaseMedFilt);

        % Plot
%         figure(7);
%         plot(timeVector(1:end-1), ressPhaseHz, 'r--'); hold on;
%         plot(timeVector(1:end-1), phaseMedFilt, 'k-'); hold on;
%         set(gca, 'xlim', [timeVector(1) timeVector(end)-1]);
%         limY = get(gca, 'ylim');
%         if strcmp(Conditions{iCondition}(end-3:end), 'Walk')
%             plot([1 length(ressPhaseHz)], [stepFreq stepFreq], 'b-');
%             if strcmp(Conditions{iCondition}(1:3), 'RAC')
%                 plot([1 length(ressPhaseHz)], [stimFreq stimFreq], 'color', [0.80,0.80,0.80]);
%                 legend({'Before moving median smoothing', 'After moving median smoothing', 'Step Frequency', 'RAC frequency'}, 'FontSize', 14);
%             else
%                 legend({'Before moving median smoothing', 'After moving median smoothing', 'Step Frequency'}, 'FontSize', 14);
%             end
%         elseif strcmp(Conditions{iCondition}(end-3:end), 'Rest')
%             plot([1 length(ressPhaseHz)], [sFreq sFreq], 'color', [0.80,0.80,0.80]);
%             if strcmp(Conditions{iCondition}(1:3), 'RAC')
%                 legend({'Before moving median smoothing', 'After moving median smoothing', 'RAC frequency'}, 'FontSize', 14);
%             else
%                 legend({'Before moving median smoothing', 'After moving median smoothing', 'RAC frequency of RAC condition '}, 'FontSize', 14);
%             end
%         end
%         xlabel({'Time (s)'}, 'FontSize', 14), ylabel({'Instantaneous Frequency (Hz)'}, 'FontSize', 14);
%         txt = (['Stability index = ' num2str(stabilityIndex)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FitBoxToText','on', 'FontSize', 14);
%         title('Instantaneous frequencies of the RESS component', 'FontSize', 16)
%         saveas(figure(7), [directoryResults '/fig_stabilityIndex.png']);
     
        % Save EEG results
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).stabilityIndex = stabilityIndex;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).phaseAngle     = ressPhaseAngle;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).phaseHz        = phaseMedFilt;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).RESS           = RESS;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).RESStime       = RESStime;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).Map            = maps(:,compMax);
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).chanLocs       = EEG.chanlocs;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).ressSNR        = snrR;
        resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).elecSNR        = snrE;
        save([pathResults 'All/' analysisType 'resultsEEG'], 'resultsEEG');

        ALLEEG = [];
        clear covarianceHi covarianceLo covarianceR covarianceS...
              data dataHi dataLo dataRaw dataS dataX elecx f FFT Hz EEG...
              lambdaIndex Lambdas lambdaSorted map2plot maps phasedMed phaseMedFilt...
              RESS ressFiltered ressHilbert ressPhaseAngle ressPhaseHz RESSsorted RESStime RESSx...
              snrE snrR temp timeVector W
        close all

    end

end