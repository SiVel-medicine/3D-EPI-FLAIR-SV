function [F0,Fn,Zn,F,info] = EPG_FLAIR(theta,ESP,T1,T2,TI,varargin)

if nargin == 0
    testrun
    F0 = [];
    Fn = [];
    Zn = [];
    F = [];
    info = [];
    return
end

opts = parse_options(varargin{:});
sequence = lower(opts.sequence);

if any(strcmp(sequence,{'epi','epi3d','3depi','3d_epi','flair_epi'}))
    [F0,Fn,Zn,F,info] = simulate_epi3d(theta,ESP,T1,T2,TI,opts);
else
    [F0,Fn,Zn,F] = simulate_tse(theta,ESP,T1,T2,TI,opts);
    info = struct('sequence','tse','echoTimes',(1:numel(F0))*ESP);
end

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
opts.nSegments = 1; 
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
        case {'nsegments','nshots','etl_per_segment'}
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

function [F0,Fn,Zn,F,info] = simulate_epi3d(theta,ESP,T1,T2,TI,opts)
alpha = abs(theta(1));
phi = angle(theta(1));
T2star = opts.T2star;
if isempty(T2star)
    T2star = T2;
end

if isempty(opts.nPE)
    if numel(theta) > 1
        nPE = numel(theta) - 1;
    else
        nPE = 64;
    end
else
    nPE = opts.nPE;
end

nPartitions = opts.nPartitions;
nSegments = max(1,round(opts.nSegments)); 
if isempty(opts.effectiveEcho)
    effectiveEcho = floor(nPE/2) + 1;
else
    effectiveEcho = opts.effectiveEcho;
end
if effectiveEcho < 1 || effectiveEcho > nPE
    error('effectiveEcho must be between 1 and nPE.');
end

echoNumber = 1:nPE;
echoTimes = echoNumber*ESP;

segBoundaries = round(linspace(0,nPE,nSegments+1)); 
segLengths = diff(segBoundaries); 

if isempty(opts.kmax)
    kmax = max(segLengths) + 2; 
else
    kmax = opts.kmax;
end

N = 3*(kmax+1);
Splus = EPG_shift_matrices(kmax);
Sminus = Splus';

kyOrder = make_k_order(nPE,opts.phaseOrdering); 
kzOrder = make_k_order(nPartitions,opts.partitionOrdering); 

F0 = zeros(1,nPE);
F = zeros(N,nPE);
segmentIndex = zeros(1,nPE); 
echoTimeInSegment = zeros(1,nPE); 

MzAfterT2Prep = zeros(1,nSegments);
MzAfterInversion = zeros(1,nSegments);
MzBeforeExcitation = zeros(1,nSegments);
MxyAfterExcitation = zeros(1,nSegments);
MzAfterExcitation = zeros(1,nSegments);

Fn = zeros(nPE,nPartitions); 

nTotalShots = nPartitions*nSegments; 
shotDurationMs = opts.T2prepTE + TI + segLengths*ESP; 
A = RF_rot(alpha,phi); 

for iz = 1:nPartitions
    isReferencePartition = (iz == 1); 
    F0part = zeros(1,nPE); 

    for seg = 1:nSegments
        nEchoSeg = segLengths(seg);
        if nEchoSeg == 0
            continue
        end

        FF = zeros(N,1);
        FF(3) = opts.zinit;

        % Dynamic FLAIR preparation in EPG states, repeated for this shot:
        % 1) optional T2 prep attenuates Mz before inversion,
        % 2) invert with a 180 degree RF pulse,
        % 3) relax during TI with T1 and T2,
        % 4) excite with the imaging flip angle.
        [FF,prepInfo] = flair_prep(FF,TI,T1,T2,kmax,opts.T2prepTE);
        FF = apply_rf(FF,A,kmax);

        if isReferencePartition
            MzAfterT2Prep(seg) = prepInfo.MzAfterT2Prep;
            MzAfterInversion(seg) = prepInfo.MzAfterInversion;
            MzBeforeExcitation(seg) = real(FF(3));
            MxyAfterExcitation(seg) = FF(1);
            MzAfterExcitation(seg) = real(FF(3));
        end

        echoRange = segBoundaries(seg)+1 : segBoundaries(seg+1); 
        for echoInSegment = 1:nEchoSeg
            FF = relax_epg(FF,0.5*ESP,T1,T2star,kmax);
            FF = Splus*FF;
            FF = relax_epg(FF,0.5*ESP,T1,T2star,kmax);
            FF = Sminus*FF;

            echo = echoRange(echoInSegment); 
            F0part(echo) = FF(1);

            if isReferencePartition
                F(:,echo) = FF;
                F0(echo) = FF(1);
                segmentIndex(echo) = seg;
                echoTimeInSegment(echo) = echoInSegment*ESP;
            end
        end
    end

    Fn(kyOrder,kzOrder(iz)) = F0part(:);
end

Zn = F(3:3:end,:);

partitionTimes = zeros(1,nPartitions);
partitionDecay = ones(1,nPartitions);

totalScanTimeMs = nPartitions*sum(shotDurationMs);

info = struct();
info.sequence = 'epi3d';
info.echoTimes = echoTimes;
info.echoNumber = echoNumber;
info.effectiveEcho = effectiveEcho;
info.effectiveTE = echoTimeInSegment(effectiveEcho); 
info.T2star = T2star;
info.T2prepTE = opts.T2prepTE;
info.nSegments = nSegments;
info.segLengths = segLengths;
info.segmentIndex = segmentIndex;
info.echoTimeInSegment = echoTimeInSegment;
info.MzAfterT2Prep = MzAfterT2Prep;
info.MzAfterInversion = MzAfterInversion;
info.MzBeforeExcitation = MzBeforeExcitation;
info.MxyAfterExcitation = MxyAfterExcitation;
info.MzAfterExcitation = MzAfterExcitation;
info.kmax = kmax;
info.FLAIRPrep = 'dynamic EPG T2 prep, inversion pulse, and TI relaxation with T1/T2, rerun independently per segment and per partition';
info.EPIReadout = 'dynamic EPG balanced dephase/rephase echo loop, restarted at each segment and partition excitation';
info.kyOrder = kyOrder;
info.kzOrder = kzOrder;
info.kyPositions = k_positions(nPE);
info.kzPositions = k_positions(nPartitions);
info.partitionTimes = partitionTimes;
info.partitionDecay = partitionDecay;
info.nTotalShots = nTotalShots; 
info.shotDurationMs = shotDurationMs; 
info.totalScanTimeMs = totalScanTimeMs; 
end

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

function FF = apply_rf(FF,A,kmax)
for k = 0:kmax
    idx = 3*k + (1:3);
    FF(idx) = A*FF(idx);
end
end

function FF = relax_epg(FF,dt,T1,T2,kmax)
E1 = exp(-dt/T1);
E2 = exp(-dt/T2);
for k = 0:kmax
    idx = 3*k + (1:3);
    FF(idx(1)) = E2*FF(idx(1));
    FF(idx(2)) = E2*FF(idx(2));
    FF(idx(3)) = E1*FF(idx(3));
end
FF(3) = FF(3) + (1-E1);
end

function order = make_k_order(n,ordering)
switch lower(ordering)
    case 'linear'
        order = 1:n;
    case 'centric'
        centre = ceil(n/2);
        offsets = 0:max(centre-1,n-centre);
        order = [];
        for ii = offsets
            candidates = [centre-ii centre+ii];
            candidates = unique(candidates,'stable');
            candidates = candidates(candidates>=1 & candidates<=n);
            order = [order candidates]; 
        end
    otherwise
        error('Unknown k-space ordering "%s".',ordering);
end
end

function k = k_positions(n)
k = -floor(n/2):ceil(n/2)-1;
end

function S = EPG_shift_matrices(kmax)
N = 3*(kmax+1);
S = sparse(N,N);

for k = 0:kmax
    fp = 3*k + 1;
    fm = 3*k + 2;
    z = 3*k + 3;

    if k < kmax
        S(fp+3,fp) = 1;
    end

    if k > 0
        S(fm-3,fm) = 1;
    else
        S(fp,fm) = 1;
    end

    S(z,z) = 1;
end
end

function E = E_diff(Ebase,diff,kmax,N)
E = sparse(zeros(N,N));
for k = 0:kmax
    idx = 3*k + (1:3);
    attenuation = exp(-diff.D*(diff.G*diff.tau*k)^2);
    E(idx,idx) = Ebase;
    E(idx(1),idx(1)) = E(idx(1),idx(1))*attenuation;
    E(idx(2),idx(2)) = E(idx(2),idx(2))*attenuation;
end
E = E(1:N,1:N);
end
