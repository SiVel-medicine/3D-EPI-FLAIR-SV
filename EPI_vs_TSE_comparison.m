
% Both readouts share the exact same flair_prep code path inside EPG_FLAIR.m ('sequence','epi3d' vs the default 'sequence','tse'), so any difference seen here is due to the readout physics itself, not a mismatched preparation:
%   - TSE: each echo is refocused by a 180 deg pulse, so signal decays with T2 (slow) across the echo train.
%   - EPI: nothing refocuses off-resonance dephasing between echoes, so signal decays with T2*

clear
clc
close all

flipAngle = deg2rad(90);   % excitation flip angle
refocusAngle = deg2rad(180); % TSE refocusing flip angle 
ESP = 8;                   % echo spacing 
ESP_epi = 0.8;              % EPI echo spacing
ETL = 64;                  % echo train length shared by both readouts for the per-shot signal comparison
T2prepTE = 50;             
TI = 3000;                 

tissue(1) = struct('name','White matter','short','WM','T1',1226,'T2',39,'T2star',32,'color',[0.000 0.278 0.671]);
tissue(2) = struct('name','Grey matter','short','GM','T1',1781,'T2',55,'T2star',42,'color',[0.800 0.157 0.157]);
tissue(3) = struct('name','CSF','short','CSF','T1',4300,'T2',1200,'T2star',450,'color',[1.000 1.000 1.000]);

%% Run TSE and EPI  for each tissue, through a single shot
tseSignal = cell(1,numel(tissue));
epiSignal = cell(1,numel(tissue));

for ii = 1:numel(tissue)
    theta = [flipAngle, repmat(refocusAngle,1,ETL-1)];
    [F0_tse,~,~,~] = EPG_FLAIR(theta,ESP,tissue(ii).T1,tissue(ii).T2,TI, ...
        'T2prepTE',T2prepTE);
    tseSignal{ii} = abs(F0_tse);

    [F0_epi,~,~,~,~] = EPG_FLAIR(flipAngle,ESP_epi,tissue(ii).T1,tissue(ii).T2,TI, ...
        'sequence','epi3d','T2star',tissue(ii).T2star,'T2prepTE',T2prepTE, ...
        'nPE',ETL,'nPartitions',1,'nSegments',1,'effectiveEcho',floor(ETL/2)+1);
    epiSignal{ii} = abs(F0_epi);
end

%% This is a figure GM/WM contrast and CSF nulling at the effective echo with a matched ETL. This comparison isolates the T2 vs T2* decay-rate difference from any readout-time difference.

effEchoIdxTSE = floor((ETL-1)/2);
effEchoIdxEPI = floor(ETL/2);

contrastTSE = tseSignal{1}(effEchoIdxTSE) - tseSignal{2}(effEchoIdxTSE); % WM-GM
contrastEPI = epiSignal{1}(effEchoIdxEPI) - epiSignal{2}(effEchoIdxEPI);
csfTSE = tseSignal{3}(effEchoIdxTSE);
csfEPI = epiSignal{3}(effEchoIdxEPI);

figure('Name','EPI vs TSE readout comparison','Position',[100 100 700 550]);
barData = [contrastTSE contrastEPI; csfTSE csfEPI];
b = bar(barData);
set(gca,'XTickLabel',{'WM-GM contrast','CSF signal'});
legend({'TSE','EPI'},'Location','best');
ylabel('|F_0| at effective echo');
title('Contrast and CSF nulling at the effective echo ');
grid on; box on;
