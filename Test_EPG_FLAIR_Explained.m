clear
clc
close all

% This is the file to run for the project.
% EPG_FLAIR.m is the function used by this test script.

flipAngle = deg2rad(90);
ESP = 0.8;                  % EPI echo spacing in ms
linesPerSegment = 64;       % new -> original code had "nPE = 64" meaning the whole (unsegmented) train length
nSegments = 3;               % new -> original code declared this but never passed it properly to EPG_flair
nPE = linesPerSegment*nSegments; % new -> old code set nPE=64 directly but now nPE is the TOTAL image resolution
nPartitions = 32;            % new -> old code set nPartitions = nSegments but rests now as a genuine independent value
readoutSamples = 64;         % points along each acquired EPI line
effectiveEcho = floor(nPE/2) + 1;
TI = 3000;                  % close to CSF null: 4300*log(2) = 2980 ms
T2prepTE = 50;               % T2 prep before inversion, in ms

tissue(1) = struct('name','White matter','short','WM','T1',1226,'T2',39,'T2star',32,'color',[0.000 0.278 0.671]);
tissue(2) = struct('name','Grey matter','short','GM','T1',1781,'T2',55,'T2star',42,'color',[0.800 0.157 0.157]);
tissue(3) = struct('name','CSF','short','CSF','T1',4300,'T2',1200,'T2star',450,'color',[1.000 1.000 1.000]);

signals = struct();  %new empty structure to store results of EPG simulations
for ii = 1:numel(tissue)
    [signals(ii).F0,signals(ii).K3D,signals(ii).Zn,signals(ii).F,signals(ii).info] = EPG_FLAIR(flipAngle,ESP,tissue(ii).T1,tissue(ii).T2,TI, 'sequence','epi3d', 'T2star',tissue(ii).T2star, 'T2prepTE',T2prepTE,'nPE',nPE, 'nPartitions',nPartitions,'nSegments',nSegments,'effectiveEcho',effectiveEcho); % new -> added 'nSegments'
end % run function 3 times for WM,GM,CSF

fprintf('Estimated total scan time: %.1f s (%.1f min) for %d shots (%d segments x %d partitions)\n', ...
    signals(1).info.totalScanTimeMs/1000, signals(1).info.totalScanTimeMs/60000, ...
    signals(1).info.nTotalShots, nSegments, nPartitions); % new -> print the timing estimate so it's not silently computed and ignored

echoNumber = signals(1).info.echoNumber;  %take echo number vector from WM simulation
effectiveTE = signals(1).info.effectiveTE; %extract effective TE (time since that shot's own excitation)
segLengths = signals(1).info.segLengths;  % new -> ky lines acquired by each of the nSegments shots (new info field)
segBoundaries = cumsum(segLengths);        % new -> used to draw segment-boundary markers in Figures 2/3/5
totalKspaceLines = nPE; % new -> it was "nPE*nSegments" but now nPE already is the total, so there is no more multiplying here
kx2d = linspace(-0.5,0.5,readoutSamples);
ky2d = linspace(-0.5,0.5,totalKspaceLines);
readoutEnvelope = exp(-18*kx2d.^2);
readoutEnvelope = readoutEnvelope/max(readoutEnvelope);
K2D = cell(1,numel(tissue));
PSF2D = cell(1,numel(tissue));
for ii = 1:numel(tissue)
    K2D{ii} = abs(signals(ii).F0(:))*readoutEnvelope; % new -> original code called fill_epi2d_kspace(F0,nSegments,readoutEnvelope), which tiled one decayed train nSegments times. I deleted the tiling function.
    psf = abs(fftshift(ifft2(ifftshift(K2D{ii}))));
    PSF2D{ii} = psf/max(psf(:));
end

wmEff = abs(signals(1).F0(effectiveEcho)); %take WM signal at effective echo
gmEff = abs(signals(2).F0(effectiveEcho));
csfEff = abs(signals(3).F0(effectiveEcho));
contrastWMGM = abs(wmEff-gmEff)/(wmEff+gmEff);  %calculate normalised contract measure between WM and GM
csfRatio = csfEff/max([wmEff gmEff]); %compare CSF signal agains whichever of GM and WM is larhger signal

fig = figure('Name','2D FLAIR EPI test','Color','w','Units','normalized','Position',[0.08 0.08 0.82 0.82]);
tabs = uitabgroup(fig);

% 1. FLAIR inversion recovery followed by the EPI readout: shows longitudinal magnetization Mz after
% inversion and the readout of the FIRST segment/shot. Every segment reruns the same FLAIR prep
% independently, so segment 1 is representative of all nSegments shots (they are identical unless
% segment-to-segment saturation is modelled).
tab = uitab(tabs,'Title','1 FLAIR + EPI');
axes('Parent',tab)
t = 0:1:TI;
hold on
hFlow = [];
flowLabels = {};
seg1Range = 1:segLengths(1); % new -> implemented because of old code having no segments, so no need to restrict to one segment's range
for ii = 1:numel(tissue)
    Mz0 = signals(ii).info.MzAfterInversion(1); % new -> old code had "signals(ii).info.MzAfterInversion" as a scalar. Now it's a per-segment vector, meaning (1) picks segment 1
    Mz = 1 + (Mz0 - 1)*exp(-t/tissue(ii).T1);
    hFlow(end+1) = plot(t,Mz,'LineWidth',2,'Color',tissue(ii).color);
    flowLabels{end+1} = sprintf('%s recovery',tissue(ii).short);
    plot(0,Mz0,'o','MarkerSize',6,'MarkerFaceColor',tissue(ii).color,'MarkerEdgeColor','w','HandleVisibility','off')
    epiTime = TI + signals(ii).info.echoTimeInSegment(seg1Range); %new -> old code used "TI + signals(ii).info.echoTimes" but it now shows segment 1's own short curve
    epiSignal = abs(signals(ii).F0(seg1Range)); % new -> old code used the full "signals(ii).F0" (whole train) but now only segment 1's echoes
    hFlow(end+1) = plot(epiTime,epiSignal,'--','LineWidth',2,'Color',tissue(ii).color);
    flowLabels{end+1} = sprintf('%s EPI (segment 1)',tissue(ii).short); % new -> old label was just "%s EPI" with no segment note
    hFlow(end+1) = plot([TI epiTime(1)],[signals(ii).info.MzBeforeExcitation(1) epiSignal(1)], ':','LineWidth',1.2,'Color',tissue(ii).color); % new -> same reason as MzAfterInversion above
    flowLabels{end+1} = sprintf('%s excitation',tissue(ii).short);
end
xline(TI,'--','excitation','LineWidth',1.4,'Color',[0.05 0.45 0.20])
xline(TI+effectiveTE,'--','k_0','LineWidth',1.4,'Color',[0.10 0.10 0.10])
yline(0,':','Color',[0.25 0.25 0.25])
yline(-1,':','no T2 prep start','Color',[0.45 0.45 0.45])
xlabel('Time after inversion (ms)')
ylabel('M_z before excitation / |F_0| during EPI')
title(sprintf('1. T2-prepared FLAIR recovery into a %d-echo segment, TE = %.1f ms (repeated identically for all %d segments)',segLengths(1),effectiveTE,nSegments)) 
legend(hFlow,flowLabels,'Location','southeast')
grid on
box on

% 2. EPI echo train shown clearly after excitation: shows signal across ALL nSegments shots
% concatenated. Each segment re-excites from equilibrium, so the signal resets and decays again
% at every segment boundary instead of decaying monotonically over the full nPE-line train.
tab = uitab(tabs,'Title','2 EPI readout');
axes('Parent',tab)
hold on
for ii = 1:numel(tissue)
    plot(echoNumber,abs(signals(ii).F0),'LineWidth',2,'Color',tissue(ii).color)
end
for ii = 1:numel(tissue) % new -> this marks the first echo of each tissue's curve
    plot(1,abs(signals(ii).F0(1)),'o','MarkerSize',5,'MarkerFaceColor',tissue(ii).color,'MarkerEdgeColor','w','HandleVisibility','off')
end
mark_segment_boundaries(segBoundaries) % new -> dotted line at every segment boundary, showing where the signal resets in terms of clarity
xline(effectiveEcho,'--','k_0','LineWidth',1.4,'Color',[0.10 0.10 0.10])
xlabel('EPI echo / ky line after excitation')
ylabel('|F_0|')
title(sprintf('2. EPI readout, %d segments of %d echoes each (dotted lines = new excitation), TE = %.1f ms',nSegments,segLengths(1),effectiveTE)) 
grid on
box on

% 3. Dynamic EPG states during the EPI readout: shows how longitudinal
% magnetisation of each tissue changes during EPI readout after FLAIR prep. The reset seen at
% each segment boundary in Figure 2 shows up here too, since Z_0 is re-inverted and re-relaxed
% independently for every segment.
tab = uitab(tabs,'Title','3 Dynamic EPG states');
axes('Parent',tab)
hold on
for ii = 1:numel(tissue)
    plot(echoNumber,real(signals(ii).F(3,:)),'LineWidth',2, ...
        'Color',tissue(ii).color)
end
mark_segment_boundaries(segBoundaries) % new -> the same principle of segment-boundary markers as Figure 2
xlabel('EPI echo number')
ylabel('Z_0 during readout')
title('3. Dynamic EPG state evolution after FLAIR preparation (dotted lines = new excitation)')
legend({tissue.short},'Location','best')
grid on
box on

% 4. Tissue signal at the k-space centre: compares the WM, GM, and CSF signal amplitudes at the effective echo, where the centre of k-space is sampled.
%sshows the resulting tissue contrast and how effectively the FLAIR preparation suppresses the CSF signal.
tab = uitab(tabs,'Title','4 Effective echo');
axes('Parent',tab)
b = bar([wmEff gmEff csfEff],'FaceColor','flat');
for ii = 1:numel(tissue)
    b.CData(ii,:) = tissue(ii).color;
end
set(gca,'XTickLabel',{tissue.short})
ylabel('|F_0| at k-space centre')
title(sprintf('4. Contrast %.3f, CSF ratio %.3f',contrastWMGM,csfRatio))
grid on
box on

% 5. 2D EPI k-space trajectory and WM/GM signal filling: shows the zigzag EPI trajectory used to
% fill k-space across the full nPE phase-encoding lines. Colour intensity comes directly from the
% per-segment simulated signal (signals(ii).F0), so each segment's block of lines now shows its
% own decay envelope restarting from a bright value instead of one train tiled nSegments times.
tab = uitab(tabs,'Title','5 2D k-space');
tl = tiledlayout(tab,1,2,'TileSpacing','compact','Padding','compact');
for ii = 1:2
    ax = nexttile(tl,ii);
    hold(ax,'on')
    values = abs(signals(ii).F0(:)); % new -> original code was "repmat(abs(signals(ii).F0(:)),nSegments,1)" which copied the same 64-value train 3 times. F0 is now already the real 192-value train, with no repmat needed.
    values = values/max(values);
    for lineIdx = 1:totalKspaceLines
        if mod(lineIdx,2) == 1
            xLine = [-0.5 0.5];
        else
            xLine = [0.5 -0.5];
        end
        yLine = ky2d(lineIdx)*[1 1];
        lineColor = 0.15 + 0.85*values(lineIdx)*tissue(ii).color;
        plot(ax,xLine,yLine,'LineWidth',1.3,'Color',lineColor)
    end
    for jj = 1:numel(segBoundaries)-1 % new -> horizontal guide is drawn at each segment boundary in k-space
        yline(ax,ky2d(segBoundaries(jj))+0.5*mean(diff(ky2d)),':','Color',[0.4 0.4 0.4],'HandleVisibility','off')
    end
    xline(ax,0,':','Color',[0.20 0.20 0.20])
    yline(ax,0,':','Color',[0.20 0.20 0.20])
    axis(ax,[-0.55 0.55 -0.55 0.55])
    xlabel(ax,'k_x readout')
    ylabel(ax,'k_y phase encode')
    title(ax,sprintf('5. %s 2D zigzag: %d segments x %d lines = %d',tissue(ii).short,nSegments,segLengths(1),totalKspaceLines)) % new -> made it so the output is real variable instead of fixed text
    grid(ax,'on')
    box(ax,'on')
end

% 6. 2D point-spread function from the filled EPI k-space:shows the 2D point-spread function produced by the simulated EPI k-space for WM and GM.
% illustrates how the k-space weighting caused by EPI signal decay affects spatial localisation and image sharpness.
tab = uitab(tabs,'Title','6 2D PSF');
tl = tiledlayout(tab,1,2,'TileSpacing','compact','Padding','compact');
for ii = 1:2
    ax = nexttile(tl,ii);
    imagesc(ax,kx2d,ky2d,PSF2D{ii})
    axis(ax,'image')
    colorbar(ax)
    colormap(ax,'parula')
    xlabel(ax,'Readout position')
    ylabel(ax,'Phase position')
    title(ax,sprintf('6. %s 2D PSF from EPI k-space',tissue(ii).short))
    grid(ax,'on')
    box(ax,'on')
end

% 7. Point-spread function from EPI T2* decay: shows PSF, EPI signal decay
% across ky causes blurring, PSF shows how much a point spreads in the
% image. Diferent tissues have different PSF due to T2star.
tab = uitab(tabs,'Title','7 PSF');
axes('Parent',tab)
hold on
for ii = 1:numel(tissue)
K = abs(signals(ii).F0);
K = K/max(K);
psf = abs(fftshift(ifft(ifftshift(K))));
psf = psf/max(psf);
plot(linspace(-0.5,0.5,numel(psf)),psf,'LineWidth',2, 'Color',tissue(ii).color)
end
xlabel('Normalised phase position')
ylabel('Normalised PSF')
title('7. EPI PSF')
legend({tissue.short},'Location','northeast')
grid on
box on

% 8. TI optimisation:evaluates how the inversion time affects WM–GM contrast and CSF suppression across the tested range.
%selected TI corresponds to the minimum simulated CSF signal, while the contrast curve shows how tissue separation changes with TI.
tab = uitab(tabs,'Title','8 TI optimisation');
axes('Parent',tab)
TI_values = 1800:25:4200;
WM = zeros(size(TI_values));
GM = zeros(size(TI_values));
CSF = zeros(size(TI_values));
contrast = zeros(size(TI_values));

for ii = 1:numel(TI_values)
    TItest = TI_values(ii);
    [F0,~,~,~] = EPG_FLAIR(flipAngle,ESP,tissue(1).T1,tissue(1).T2,TItest, ...
        'sequence','epi3d','T2star',tissue(1).T2star,'T2prepTE',T2prepTE, ...
        'nPE',nPE,'nPartitions',nPartitions,'nSegments',nSegments,'effectiveEcho',effectiveEcho); % new -> added 'nSegments'
    WM(ii) = abs(F0(effectiveEcho));

    [F0,~,~,~] = EPG_FLAIR(flipAngle,ESP,tissue(2).T1,tissue(2).T2,TItest, ...
        'sequence','epi3d','T2star',tissue(2).T2star,'T2prepTE',T2prepTE, ...
        'nPE',nPE,'nPartitions',nPartitions,'nSegments',nSegments,'effectiveEcho',effectiveEcho); % same reason as above
    GM(ii) = abs(F0(effectiveEcho));

    [F0,~,~,~] = EPG_FLAIR(flipAngle,ESP,tissue(3).T1,tissue(3).T2,TItest, ...
        'sequence','epi3d','T2star',tissue(3).T2star,'T2prepTE',T2prepTE, ...
        'nPE',nPE,'nPartitions',nPartitions,'nSegments',nSegments,'effectiveEcho',effectiveEcho); % new -> same reason as above
    CSF(ii) = abs(F0(effectiveEcho));
    contrast(ii) = abs(WM(ii)-GM(ii))/(WM(ii)+GM(ii));
end

[~,csfMinIdx] = min(CSF);
bestNullTI = TI_values(csfMinIdx);

purple = [0.45 0.10 0.60];
black = [0.05 0.05 0.05];
green = [0.05 0.45 0.20];

yyaxis left
hContrast = plot(TI_values,contrast,'LineWidth',2.2,'Color',purple);
ylabel('WM-GM contrast')
ax = gca;
ax.YAxis(1).Color = purple;

yyaxis right
hCSF = plot(TI_values,CSF,'LineWidth',2.2,'Color',black);
ylabel('CSF signal')
ax.YAxis(2).Color = black;

hNull = xline(bestNullTI,'--','minimum CSF','LineWidth',1.4,'Color',green);
xlabel('TI (ms)')
title(sprintf('8. TI optimisation, best TI = %.0f ms',bestNullTI))
legend([hContrast hCSF hNull],{'Purple: WM-GM contrast', ...
    'Black: CSF signal','Green: best CSF null'},'Location','best')
grid on
box on


% new -> this file used to end with a function fill_epi2d_kspace(echoTrain,nSegments,readoutEnvelope) that copy-pasted an echo train nSegments times which may lead to a false segmented k-space fill. 
% It is no longer needed, with signals(ii).F0 now being the real, per-segment train.

function mark_segment_boundaries(segBoundaries) %new -> this draws a dotted vertical guide at every segment boundary on the current axes.
% Used by Figures 2 and 3 primarily
for ii = 1:numel(segBoundaries)-1
    xline(segBoundaries(ii)+0.5,':','Color',[0.5 0.5 0.5],'HandleVisibility','off')
end
end
