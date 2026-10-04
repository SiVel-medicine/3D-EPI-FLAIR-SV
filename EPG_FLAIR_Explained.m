function [F0,Fn,Zn,F,info] = EPG_FLAIR(theta,ESP,T1,T2,TI,varargin)

% Default behaviour keeps the original TSE EPG model:
% [F0,Fn,Zn,F] = EPG_FLAIR(theta,ESP,T1,T2,TI)
% 3D FLAIR EPI mode models an inversion-prepared excitation followed by a
% gradient echo EPI readout:
% [F0,Fn,Zn,F,info] = EPG_FLAIR(deg2rad(90),ESP,T1,T2,TI,
% 'sequence','epi3d','T2star',45,'nPE',64,'nPartitions',32)
% instead of theta=[] use theta=deg2rad(90) as EPI only has one RF pulse
% during imaging as opposed tp TSE
% For EPI, F0 is the echo train sampled along ky, Fn is the 3D ky-kz
% weighting matrix, Zn is the longitudinal state after excitation, and F
% stores the compact EPG state [F+0; F-0; Z0] over the readout.

if nargin == 0 % if number of input arguments=0 then run testrun.m
testrun
F0 = [];
Fn = [];
Zn = [];
F = [];
info = [];
return %return empty outputs so MATLAB doesnt complain outputs werent assigned
end

opts = parse_options(varargin{:}); %varargin is optional user input so sequence,epi etc are stored inside varargin
sequence = lower(opts.sequence); %converts EPI3D into epi3D so comparisons are easier

if any(strcmp(sequence,{'epi','epi3d','3depi','3d_epi','flair_epi'}))
[F0,Fn,Zn,F,info] = simulate_epi3d(theta,ESP,T1,T2,TI,opts);
else
[F0,Fn,Zn,F] = simulate_tse(theta,ESP,T1,T2,TI,opts);
info = struct('sequence','tse','echoTimes',(1:numel(F0))*ESP);
end
%if user request EPI do it, otherwise proceed with TSE
end

function opts = parse_options(varargin)
if mod(numel(varargin),2) ~= 0
error('Optional inputs must be name/value pairs.');
end

opts = struct();
opts.sequence = 'tse';
opts.kmax = [];
opts.diff = [];
opts.zinit = 1;
opts.T2star = [];
opts.T2prepTE = 0;
opts.nPE = [];
opts.nPartitions = 32;
opts.nSegments = 1; % new -> number of independent shots the nPE ky lines are split across
opts.effectiveEcho = [];
opts.phaseOrdering = 'linear';
opts.partitionOrdering = 'linear';
opts.partitionSignalDecay = false;
opts.partitionTR = [];

for ii = 1:2:numel(varargin)
name = lower(varargin{ii});
value = varargin{ii+1};
switch name
case 'sequence'
opts.sequence = value;
case 'kmax'
opts.kmax = value;
case 'diff'
opts.diff = value;
case 'zinit'
opts.zinit = value;
case {'t2star','t2s','t2prime'}
opts.T2star = value;
case {'t2prep','t2prepte','t2preparation'}
opts.T2prepTE = value;
case {'npe','etl','echotrainlength'}
opts.nPE = value;
case {'npartitions','nz','nslices'}
opts.nPartitions = value;
case {'nsegments','nshots','etl_per_segment'} % new -> allows to call nsegments
opts.nSegments = value;
case {'effectiveecho','centreecho','centerecho'}
opts.effectiveEcho = value;
case {'phaseordering','kyordering'}
opts.phaseOrdering = value;
case {'partitionordering','kzordering'}
opts.partitionOrdering = value;
case 'partitionsignaldecay'
opts.partitionSignalDecay = logical(value);
case {'partitiontr','shottr'}
opts.partitionTR = value;
otherwise
error('Unknown EPG_FLAIR option "%s".',varargin{ii});
end
end
end
% New function simulate_epi3d
%receives: theta=RF flip angles, echo spacing, T1, T2, TI, additional seq
%parameters(previous varargin)
%returns: F0=signal at every EPI echo, Fn=signal mapped into 3D k-space,
%F=all EPG states throughout readout, info=useful sequence info
function [F0,Fn,Zn,F,info] = simulate_epi3d(theta,ESP,T1,T2,TI,opts)
alpha = abs(theta(1)); %extracts flip angle of excitation pulse
phi = angle(theta(1)); %extracts phase of RF pulse, usually ignired but EPG needs it to determine which coherence pathways are created, sign of transverse magnetisation, phase cycling
T2star = opts.T2star; %T2star contains true T2+magnetic field inhomogeneity+susceptibility+off-resonance, T2star<T2 hence signal decays faster
if isempty(T2star)
T2star = T2;
end

if isempty(opts.nPE) %is number of phase encoding lines specified?
if numel(theta) > 1 %checking if user accidentally passed TSE flip angle train
nPE = numel(theta) - 1; %number of phase encoding lines
else
nPE = 64; % if there is only 1 RF pulse assume 64 ky lines; 64 is random, real scanner can be more 
end
else
nPE = opts.nPE; %can supply any value
end
%unique to 3D EPI, ky x kx x kz so need kz=partition encoding
%instead of exciting one slice sequence excites entire slab and separates
%it into many partitions
nPartitions = opts.nPartitions;
nSegments = max(1,round(opts.nSegments)); % new -> valeriia's code had no concept of segments at all. Allows at least 1 shot with nSegments>1 spliting nPE ky lines across excitations
if isempty(opts.effectiveEcho) %was which echo acquires centre of k space specified?
effectiveEcho = floor(nPE/2) + 1; %if not specified suppose 64 echoes/2+1
else
effectiveEcho = opts.effectiveEcho; %if specified just use that
end
if effectiveEcho < 1 || effectiveEcho > nPE %make sure effctive echo exists
error('effectiveEcho must be between 1 and nPE.');
end

echoNumber = 1:nPE;
echoTimes = echoNumber*ESP; %acquisition times of every gradient echo, counted continuously across all segments

% Long echo train loses signal to decay by the end, since every echo comes from a single excitation. 
% Real segmented EPI splits the lines across several shots (each with its own excitation) so each shot's train is shorter and decays less.

segBoundaries = round(linspace(0,nPE,nSegments+1)); % new ->  ky line index where each segment starts/ends
segLengths = diff(segBoundaries); % new -> number of echoes acquired by each segment

if isempty(opts.kmax)
kmax = max(segLengths) + 2; % new -> sized for the longest SINGLE shot, with every shot resets FF fresh and only ever builds up to its own echo count of coherence orders
else
kmax = opts.kmax;%if kmax specifies then just use it
end
%number of EPG variables
N = 3*(kmax+1); %every order contains F+,F-,Z and +1 accounts for F0
Splus = EPG_shift_matrices(kmax); %matrix performing the gradient pulse shifting magnetization into higher coherence order (moved from a later point)
Sminus = Splus'; % ' is matrix transpose , shifts coherence orders in opposite directions (moved from a later point)

kyOrder = make_k_order(nPE,opts.phaseOrdering); % new -> moved up from below, needed inside the partition loop to work i believe
kzOrder = make_k_order(nPartitions,opts.partitionOrdering); % same reason

F0 = zeros(1,nPE);%creates vector staring the observable MRI signal for every echo (moved from a later point)
F = zeros(N,nPE); %creates matrix that stores the complete EPG state after every echo (moved from a later point)
segmentIndex = zeros(1,nPE); % new -> which segment is acquired in each echo
echoTimeInSegment = zeros(1,nPE); % new -> time since this segment's own excitation, not since echo 1

MzAfterT2Prep = zeros(1,nSegments); % new -> was a scalar with one prep and a whole train, and is now a per-segment 
MzAfterInversion = zeros(1,nSegments); % same reason
MzBeforeExcitation = zeros(1,nSegments); % same reason
MxyAfterExcitation = zeros(1,nSegments); % same reason
MzAfterExcitation = zeros(1,nSegments); % same reason

Fn = zeros(nPE,nPartitions); % new -> moved up so the partition loop below can fill it directly

nTotalShots = nPartitions*nSegments; % total number of independent excitations in the whole acquisition
shotDurationMs = opts.T2prepTE + TI + segLengths*ESP; % per-segment shot duration (T2prep+TI+readout), vectorized once rather than recomputed inside the loop below
A = RF_rot(alpha,phi); %excitation rotation matrix does not depend on seg/iz, so build it once instead of every shot

% Each partition obtain its own fresh FLAIR prep & excitation, instead of reusing one ky train for every partition. 
% Within a partition, the ky lines are still split across the segment shots. 
% Once that partition is done, the next one starts from a new excitation

for iz = 1:nPartitions % new -> outer partition loop
isReferencePartition = (iz == 1); %single source of truth for the "only partition 1 writes the 2D reference outputs" gate used below
F0part = zeros(1,nPE); % new -> this partition's own ky signal train

for seg = 1:nSegments %run one full FLAIR prep + EPI readout per segment/shot, with a change -> original code had no outer loop, prep+excite+decay ran ONCE for the whole train
nEchoSeg = segLengths(seg);
if nEchoSeg == 0 % new -> implemented as a guard for a segment with 0 lines, in cases if nSegments > nPE
continue
end

FF = zeros(N,1); %create column vector % new -> edited the original code that created FF once, outside any loop
FF(3) = opts.zinit; %set equilibrium longitudinal magnetisation for this shot

% Dynamic FLAIR preparation in EPG states, repeated independently per segment:
% 1. optional T2 prep attenuates Mz before inversion.
% 2. 180° inversion pulse → creates the inverted longitudinal magnetisation.
% 3. T2/T1 relaxation during TI → magnetisation evolves during the inversion time. The transverse components decay according to T2, while longitudinal recovery occurs according to T1.
% 4. Imaging excitation flip angle → converts the remaining longitudinal magnetisation into transverse magnetisation.
[FF,prepInfo] = flair_prep(FF,TI,T1,T2,kmax,opts.T2prepTE); % new -> code is nowa allowed to be called once PER segment PER partition
FF = apply_rf(FF,A,kmax); %applying RF pulse to every EPG coherence state % new -> original code was excited ONCE total but excites fresh magnetization every segment/partition

if isReferencePartition % new -> describing partition 1 as the reference readout, only partition 1 writes these
MzAfterT2Prep(seg) = prepInfo.MzAfterT2Prep; % new -> indexed by seg
MzAfterInversion(seg) = prepInfo.MzAfterInversion; % same reason
MzBeforeExcitation(seg) = real(FF(3)); %stores longitudinal magnetization immediatelly before excitation % same reason
MxyAfterExcitation(seg) = FF(1); % indexed by seg
MzAfterExcitation(seg) = real(FF(3)); % same reason
end

echoRange = segBoundaries(seg)+1 : segBoundaries(seg+1); %global echo indices belonging to this segment 
%  new -> This loop now runs only over the echoes in the current segment, not all of them. 
% The old code looped over every echo in one long train, so the signal decayed the whole way through.
FF = relax_epg(FF,0.5*ESP,T1,T2star,kmax); %now use T2*, not T2, because gradient echoes are sensitive to field inhomogeneities.
FF = Splus*FF; %gradient shifts coherenec order
FF = relax_epg(FF,0.5*ESP,T1,T2star,kmax); %magnetixation continues decaying while travelling towards echo
FF = Sminus*FF; %opposite gradient applied rephasing spins forming gradient echo

echo = echoRange(echoInSegment); % new -> maps this segment's local echo count to its ky line index
F0part(echo) = FF(1); % new -> recorded into THIS partition's own array, not a shared F0

if isReferencePartition % new -> reference F/F0/segmentIndex/echoTimeInSegment outputs only recorded from partition 1
F(:,echo) = FF; %complete EPG state vector stored in corresponding column of F
F0(echo) = FF(1);%record signal amplitude of current echo
segmentIndex(echo) = seg; % new
echoTimeInSegment(echo) = echoInSegment*ESP; % new
end
end
end

Fn(kyOrder,kzOrder(iz)) = F0part(:); % new -> this partition's OWN independently-simulated signal, not F0 scaled by a post-hoc factor
end

Zn = F(3:3:end,:); %extracts all longitudinal states for every echo

%new -> The old after-the-fact partition decay scaling is no longer needed, since every partition now has its own excitation and decay. 
% new -> The opts.partitionSignalDecay and opts.partitionTR options have no effect.
partitionTimes = zeros(1,nPartitions);
partitionDecay = ones(1,nPartitions);

% Total scan time estimate: each partition repeats the same nSegments
% shot durations. This is arithmetic on quantities already computed
% above, not a new simulation.
totalScanTimeMs = nPartitions*sum(shotDurationMs);

info = struct();
info.sequence = 'epi3d';
info.echoTimes = echoTimes;
info.echoNumber = echoNumber;
info.effectiveEcho = effectiveEcho;
info.effectiveTE = echoTimeInSegment(effectiveEcho); % new -> was echoTimes(effectiveEcho), now time since that echo's own segment excitation
info.T2star = T2star;
info.T2prepTE = opts.T2prepTE;
info.nSegments = nSegments; % new
info.segLengths = segLengths; % new
info.segmentIndex = segmentIndex; % new
info.echoTimeInSegment = echoTimeInSegment; % new
info.MzAfterT2Prep = MzAfterT2Prep; % new ->  prepInfo.MzAfterT2Prep (scalar), is now per-segment vector
info.MzAfterInversion = MzAfterInversion; % new same reason
info.MzBeforeExcitation = MzBeforeExcitation; % new same reason
info.MxyAfterExcitation = MxyAfterExcitation; % new same reason
info.MzAfterExcitation = MzAfterExcitation; % new same reason
info.kmax = kmax;
info.FLAIRPrep = 'dynamic EPG T2 prep, inversion pulse, and TI relaxation with T1/T2, rerun independently per segment/partition';
info.EPIReadout = 'dynamic EPG balanced dephase/rephase echo loop, restarted at each segment & partition excitation';
info.kyOrder = kyOrder;
info.kzOrder = kzOrder;
info.kyPositions = k_positions(nPE);
info.kzPositions = k_positions(nPartitions);
info.partitionTimes = partitionTimes;
info.partitionDecay = partitionDecay;
info.nTotalShots = nTotalShots; % new
info.shotDurationMs = shotDurationMs; % new 
info.totalScanTimeMs = totalScanTimeMs; % new 
end
%%%tse
function [F0,Fn,Zn,F] = simulate_tse(theta,ESP,T1,T2,TI,opts)
np = length(theta);
if isempty(opts.kmax)
kmax = 2*(np - 1);
else
kmax = opts.kmax;
end

if isinf(kmax)
allpathways = true;
kmax = 2*(np - 1);
else
allpathways = false;
end

if allpathways
kmax_per_pulse = 2*(1:np);
kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
else
kmax_per_pulse = 2*[1:ceil(np/2) (floor(np/2)):-1:1]+1;
kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
if max(kmax_per_pulse)<kmax
kmax = max(kmax_per_pulse);
end
end

N = 3*(kmax+1);
alpha = abs(theta);
phi = angle(theta);
phi(2:end) = phi(2:end) + pi/2;
S = sparse(EPG_shift_matrices(kmax));
E1 = exp(-0.5*ESP/T1);
E2 = exp(-0.5*ESP/T2);
b = zeros([N 1]);
b(3) = 1-E1;

if ~isempty(opts.diff)
E = E_diff(diag([E2 E2 E1]),opts.diff,kmax,N);
else
E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N);
end

SE = sparse(S*E);
T = sparse(zeros(N,N));
i1 = [];
for ii = 1:3
i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3)));
end

F = zeros([N np-1]);
FF = zeros([N 1]);
FF(3) = opts.zinit;

if TI > 0 || opts.T2prepTE > 0
    [FF,~] = flair_prep(FF,TI,T1,T2,kmax,opts.T2prepTE);
end

A = RF_rot(alpha(1),phi(1));
FF(1:3) = A*FF(1:3);
kidx = 1:6;
FF(kidx) = SE(kidx,kidx)*FF(kidx)+b(kidx);

for jj = 2:np
A = RF_rot(alpha(jj),phi(jj));
build_T(A);
kidx = 1:3*kmax_per_pulse(jj);
FF(kidx) = T(kidx,kidx)*FF(kidx);
F(kidx,jj-1) = SE(kidx,kidx)*FF(kidx)+b(kidx);
F(1,jj-1) = conj(F(1,jj-1));

if jj == np
break
end
FF(kidx) = SE(kidx,kidx)*F(kidx,jj-1)+b(kidx);
FF(1) = conj(FF(1));
end

F0 = F(1,:);
idx = [fliplr(5:3:size(F,1)) 1 4:3:size(F,1)];
kvals = -kmax:kmax;
idx(1:2) = [];
kvals(1:2) = [];
Fn = F(idx,:);
Fn(kvals<0,:) = conj(Fn(kvals<0,:));
Zn = F(3:3:end,:);

function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2 = 1:9
            T(i1(i2):ksft:end) = AA(i2);
        end
    end
end

function [FF,prepInfo] = flair_prep(FF,TI,T1,T2,kmax,T2prepTE)
if nargin < 6 || isempty(T2prepTE)
    T2prepTE = 0;
end
if T2prepTE < 0
    error('T2prepTE must be zero or positive.');
end

if T2prepTE > 0
    T2prepDecay = exp(-T2prepTE/T2);
    for k = 0:kmax
        idx = 3*k + (1:3);
        FF(idx(3)) = T2prepDecay*FF(idx(3));
    end
end
prepInfo.MzAfterT2Prep = real(FF(3));

Ainv = RF_rot(pi,0);
FF = apply_rf(FF,Ainv,kmax);
prepInfo.MzAfterInversion = real(FF(3));

if TI > 0
    FF = relax_epg(FF,TI,T1,T2,kmax);
end
end

function Tap = RF_rot(a,p)
Tap = zeros([3 3]);
Tap(1) = cos(a/2).^2;
Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
Tap(4) = conj(Tap(2));
Tap(5) = Tap(1);
Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
Tap(7) = -1i*exp(1i*p)*sin(a);
Tap(8) = 1i*exp(-1i*p)*sin(a);
Tap(9) = cos(a);
end
%%%
function FF = apply_rf(FF,A,kmax) % FF=complete EPG state vector, A=3x3 RF rotation matrix, kmax=highest coherence order being simulated
for k = 0:kmax
idx = 3*k + (1:3); %determines which 3 positions in FF correspond to current coherence order
FF(idx) = A*FF(idx); %multiple the 3 states by rotation matrix
end
end

function FF = relax_epg(FF,dt,T1,T2,kmax) %what happens to EPG states during dt where there is no RF pulse
E1 = exp(-dt/T1); %calculates T1 relaxation factor
E2 = exp(-dt/T2);
for k = 0:kmax %goes through every EPG coherence order
idx = 3*k + (1:3); %identifies F+,F-,Z
FF(idx(1)) = E2*FF(idx(1));
FF(idx(2)) = E2*FF(idx(2));
FF(idx(3)) = E1*FF(idx(3));
end
FF(3) = FF(3) + (1-E1); %T1 recovery
end
%k-space acquisition
function order = make_k_order(n,ordering) %determines order in which phase-encoding lines are filled
switch lower(ordering) %convert requested ordering to lowercase and determine which acquisition startegy to use
case 'linear'
order = 1:n;
    case 'centric' %centre of k-space is acquired first or close to first
centre = ceil(n/2);
offsets = 0:max(centre-1,n-centre);
order = [];
for ii = offsets
candidates = [centre-ii centre+ii];
candidates = unique(candidates,'stable'); %remove duplicates
candidates = candidates(candidates>=1 & candidates<=n); %remove invalid indices
order = [order candidates]; %add new candidates to acquisition order
end
otherwise
error('Unknown k-space ordering "%s".',ordering);
end
end

function k = k_positions(n) %creates actual k-space coordinate values
k = -floor(n/2):ceil(n/2)-1;
end

function S = EPG_shift_matrices(kmax) %creates matrix that represents a gradient-induced coherence shift
N = 3*(kmax+1); %matrix size
S = sparse(N,N);

for k = 0:kmax %loop through coherence orders
fp = 3*k + 1; %where is F+ stored
fm = 3*k + 2;
z = 3*k + 3;

if k < kmax %positive coherence shift
S(fp+3,fp) = 1;
end

if k > 0 %negative coherence shift
S(fm-3,fm) = 1;
else
S(fp,fm) = 1; %special k=0 case
end

S(z,z) = 1; %gradient doesnt change coherence order of longitudinal magnetisation
end
end

function E = E_diff(Ebase,diff,kmax,N) %creates modified relaxation matrix that includes diffusion attenuation
E = sparse(zeros(N,N)); %create matrix
for k = 0:kmax %go through every EPG order
idx = 3*k + (1:3); %idx= [F+ F- Z]
attenuation = exp(-diff.D*(diff.G*diff.tau*k)^2); %S=S0e^-bD where b=diffusion weighting factor, D=diffusion coefficient
E(idx,idx) = Ebase; %copies normal relaxation matrix into current EPG block
E(idx(1),idx(1)) = E(idx(1),idx(1))*attenuation; 
E(idx(2),idx(2)) = E(idx(2),idx(2))*attenuation;
end
E = E(1:N,1:N); %final matrix sizing
end
