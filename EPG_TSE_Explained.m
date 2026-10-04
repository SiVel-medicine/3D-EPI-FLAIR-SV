clear
clc
addpath(genpath('EPGX_functions'))

% simulates TSE by tracking how magnetisation moves through different states
% after RF excitation pulse, 180 deg refocusing pulses, relaxation,
% gradients and predicts the signal at each echo.
% F0=MRI signal intensity at every echo
% Fn=transverse magnetisation states
% Zn=longitudinal magnetisation states
% F=complete EPG state matrix
% theta=RF pulse train
% ESP=echo spacing, eg 10 ms
% T1=longitudinal relaxation
% T2=transverse relaxation
% varargin=optional inputs later, eg kmax
% kmax=how many EPG states to calculate, eg 'kmax', 20 keep F0-F20, Z0-Z20

function [F0,Fn,Zn,F] = EPG_TSE(theta,ESP,T1,T2,varargin)

% this section checks if I gave extra parameters
for ii=1:length(varargin)
    
    if strcmpi(varargin{ii},'kmax')  %STRing CoMPare Ignore case which compares 2 strings and is true enters the if statement
        kmax = varargin{ii+1};
    end
% normally EPG assumes no diffusion, but in reality there is diffusion
% structure: gradient amplitude( diff.G), gradient duration( diff.tau),
% diffusion coefficient( diff.D), so if i want diffusion i create all these
% variables and pass the structure into function
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
% zinit tells simulator what the initial longirudinal magnetisation should
% be before the first excitation pulse, by default it assumes fully relaxed
% tissue 1, but can use another value, particularly important for FLAIR as
% after the inversion pulse Mz=-1 instead of 1, then during TI it slowly
% recovers for when the excitation pulse arrives
    if strcmpi(varargin{ii},'zinit')
        zinit=varargin{ii+1};
    end
% adding TI for FLAIR sequence
    if strcmpi(varargin{ii},'TI')
    TI = varargin{ii+1};
end
    
end



np = length(theta); % counts RF pulses
if ~exist('kmax','var') % if the variable kmax doesn't exist
    kmax = 2*(np - 1);  % (np-1) because first pulse is excitation, so only np-1 are refocusing pulses which increase EPG order
end % 2*(np-1) as every refocusing pulse generates new coherence pathways so after n refocusing pulses largest order is 2*n

if isinf(kmax) % is kmax infinite?
    allpathways = true; % putting inf in EPG_TSE after 'kmax' is used top calculate evry possible coherence pathway and useful for checking that no signal pathway is discarded
    kmax = 2*(np - 1); % interpret infinity as use the largest physically possible number of pathways as we need to create matrices and cant create infinity matrices
else
    allpathways = false; % means we are allowed to optimise and later the program will ignoire coherence pathways that are impossible to populate to speed up simulation
end




if allpathways
    kmax_per_pulse = 2*(1:np); % after each RF pulse allow more coherent pathways
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax; % only select entries where condition kmax_per_pulse>kmax is true so maximum order always <= kmax
else
    kmax_per_pulse = 2*[1:ceil(np/2) (floor(np/2)):-1:1]+1; % ceil=round up, floor=round down, middle -1 is the spet size, so decreses by 1 until reaches 1, modelling exactly the TSE echo train: at first few coherence pathways exist, more RF pulses applied and more pathways created, near middle of the sequence number is the largest and towards the end the pathways have relaxed
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax; % +1 ensures zero-order coherence is included in allocation
     
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end


N=3*(kmax+1); % +1 accounts for zero-order coherence, every coherence order contains three independent magnetisation components( F+k, F-k, Zk)
alpha = abs(theta); % how much do i rotate magnetisation?
phi=angle(theta); % angle() returns phase of a complex number
% add CPMG phase( technique that measures T2 relaxation times, initil 90
% followed by train of 180 refocusing pulses 
phi(2:end) = phi(2:end) + pi/2; % phi(2:end) means modify only refocusing pulses, because excitation pulse establishes the initial transverse magnetisation and refocusing pulses control the scho formation
    

S = EPG_shift_matrices(kmax); % create mathematical operation that represents gradient induced dephasing
S = sparse(S); % sparse matrix stores only non-zero values
% Relaxation over half the echo spacing, EPG divides echo spacing into 2
% pieces: RF pulse-half ESP-echo-half ESP-next RF pulse
E1 = exp(-0.5*ESP/T1); % longitudinal recovery factor
E2 = exp(-0.5*ESP/T2); % transverse decay 
E = diag([E2 E2 E1]); % creates relaxation matrix
b = zeros([N 1]); % creates a column vector 
b(3) = 1-E1; % index 3 is Z0 because only zero-order longitudinal state represents uniform magnetisation


if exist('diff','var') % did user provide diffusion parameters?
    E = E_diff(E,diff,kmax,N); % takes relaxation matrix E and modifies it with diffusion
else
    % If no diffusion, E is the same for all EPG orders
    E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N); % repmat repeats calculation and , makes it a column and spdiags puts values on diagonal creating a specific matrix
end
    
SE=S*E; % first E decays magnetization then S shifts coherence so SE combines relaxation and gradient dephasing
SE=sparse(SE); % only store non zero values
T = zeros(N,N); % creates empty matrix
T = sparse(T); % most elements are zero so convert to sparse

i1 = []; % cretaes empty vector
for ii=1:3 % loops 3 times as EPG order has 3 states 
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3))); % where are the first 3by3 blocks located?; sub2ind converts (row,column) into (single matrix index)
end


F = zeros([N np-1]); % F stores snapshots of magnetisation at every echo
Mz_init = 1; % default fully relaxed magnetization
if exist('TI','var') % calculate longitudinal magnetisation after inversion recovery(FLAIR)
    Mz_init = 1 - 2*exp(-TI/T1);
end
FF = zeros([N 1]); % FF is current magnetisation state changing continuosly
if ~exist('zinit','var') % is there an initial longitudinal magnetisation?
    FF(3)=zinit; % FF(3) means Z0(longitudinal magnetization immediately before excitation) 
else
    FF(3)=Mz_init; % instead of assuming Mz=1, start from whatever value I provide
end
    
% now simulation actually applies first RF pulse
A = RF_rot(alpha(1),phi(1)); %creates RF rotation matrix for first pulse
kidx = 1:3; % defines which entries of the state vector will be affected
FF(kidx) = A*FF(kidx); % RF rotation matrix A rotates this magnetisation


%now simulate the dephase gradient & evolution, half the readout.
%Immediately after excitation spins are approximately in phase
%corresponding to F0, the observable signal. Now when gradient is applied
%magnetic field is no longer identical everywhere and spins have dephased,
%moving magnetisation from F0 towards F1, etc. During evolution 3 things
%occur simltaneously: T2 decay, T1 recovery, gradient dephasing changing
%coherence orders. Matrix SE cantains all 3 together. The simulator divides
%each echo spacing into 2 equal parts allowing echo to oocur exactly in the
%middle, that's why half the readout.
kidx=1:6; % six entries: F0, F0*, Z0, F1, F-1*, Z1. Include F1 and Z1 as immediately after excitation first gradient cretaes first-order coherence
FF(kidx) = SE(kidx,kidx)*FF(kidx)+b(kidx); % FF(kidx)=replace current magnetisation, SE(kidx,kidx)=contains T2 decay, T1 relaxation, gradient shifts, kidx,kidx means take only upper left 6x6 of the matrix
%*FF(kidx) means allow the magnetisation to evolve, +b(kidx) adds
%longtidunal regrowth. Physically this line means take the current
%magnetisation, let it undergo T1 recovery, T2 decay and gradient
%dephasing for half an echo spacing, then store the updated magnetisation.

%now simulate the refocusing pulses
for jj=2:np % pulse starts at 2
    A = RF_rot(alpha(jj),phi(jj)); % every iteration the code builds RF rotation matrix for the particular pulse 
    build_T(A); % replicate A to make large T matrix
    kidx = 1:3*kmax_per_pulse(jj); % how many coherence orders could physically exist after this particular pulse
    FF(kidx)=T(kidx,kidx)*FF(kidx); % application of refocsing RF pulse 
    F(kidx,jj-1) = SE(kidx,kidx)*FF(kidx)+b(kidx); % F=snapshots at echo times, here again T2 decay, T   recovery, gradient evolution, coherence shifting are stored in SE
    F(1,jj-1)=conj(F(1,jj-1)); % 1 is first observable signal so F0, F0 comes from F-1 so have to take conjugate
    
    if jj==np
        break % so if there is no next RF pulse after final echo no need to evolve further to save computational time
    end
    
    FF(kidx) = SE(kidx,kidx)*F(kidx,jj-1)+b(kidx); % *F instead of FF as the echo has just been recorded so need to start evolving from echo to next RF pulse
    FF(1)=conj(FF(1)); % F0 comes from F-1 so conjugate 
    
end

%produce the quantity MRI scanner actually easures as entire TSE sequence
%has been simulated, RF pulse applied, echo has been generated, EPG state
%has evolved
F0 = F(1,:)*1i; % echo train, multiply by 1i as a phase convention
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)]; % this line creates correct ordering for plotting
kvals = -kmax:kmax; % cretaes labels for each coherent order
idx(1:2)=[];
kvals(1:2)=[]; % remove lowest 2 negative states since they are never populated


Fn = F(idx,:); % creates the final transverse-state matrix
Fn(kvals<0,:)=conj(Fn(kvals<0,:)); % negative coherence pathways must be conjugated
Zn = F(3:3:end,:); % extracts all longitudinal states without reordering


% defining 2 helper functions RF_rot and build_T. RF_rot computes 3x3
% rotation matrix for a pulse with flip angle a and phase p. build_T takes
% the small 3x3 matrix A from RF_rot and copies it along the diagonal of
% the large NxN matrix T, so that same RF rotation is applied independently
% to every coherence order
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

    function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end