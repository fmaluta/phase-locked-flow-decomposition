function out = flow_decomposition(DATA, META, OPT)
%FLOW_DECOMPOSITION Phase-locked flow decomposition and energy triad analysis.
%
%   out = FLOW_DECOMPOSITION(DATA, META, OPT) performs a phase-locked
%   decomposition of a time-resolved velocity field on an x–z plane and
%   computes a plane-averaged energy budget separating mean, periodic
%   (phase-locked), and residual components.
%
%   The method identifies coherent harmonics of the impeller shaft frequency
%   from a spectral analysis of a plane-averaged kinetic-energy proxy, and
%   reconstructs the associated velocity contribution via linear regression
%   onto sine/cosine bases at the selected frequencies.
%
%   The routine is designed for CFD datasets of stirred tanks or single-use
%   bioreactors with a bottom-mounted impeller, but is otherwise geometry-
%   agnostic.
%
% -------------------------------------------------------------------------
% INPUTS
%
%   DATA   struct with fields:
%       t   [T x 1]           time vector (uniformly sampled)
%       X   [Nz x Nx]         x-coordinate grid of the plane
%       Z   [Nz x Nx]         z-coordinate grid of the plane
%       U   [Nz x Nx x T]     velocity component in x-direction
%       V   [Nz x Nx x T]     velocity component in y-direction
%       W   [Nz x Nx x T]     velocity component in z-direction
%       S   [Nz x Nx] or
%           [Nz x Nx x T]     optional mask (1 = valid sample, 0 = excluded)
%
%   META   struct with fields:
%       rpm        impeller rotational speed [rpm]
%       nblades    number of impeller blades
%       name       case name (used in figure titles)
%
%   OPT    struct with optional fields (all have defaults):
%
%     OPT.compute
%       .triad_budget   compute and plot energy triad (default: true)
%       .save_figs      save figures to disk (default: true)
%       .show           display figures on screen (default: true)
%
%     OPT.harmonics
%       .Nmax           maximum harmonic index n (n*f_shaft <= Nyquist)
%       .min_frac       minimum normalized contrast to retain a harmonic
%       .keep_shaft     retain shaft frequency (default: true)
%       .guard_Hz       exclusion half-band around shaft if keep_shaft=false
%       .force_lines    optional additional frequencies (Hz or integer
%                       multiples of shaft frequency)
%
%     OPT.welch
%       .Tw             Welch window length [s]
%       .ovlp           window overlap fraction
%       .nfft           FFT length (default: nextpow2)
%       .minL           minimum window length in samples
%
%     OPT.psd
%       .highpass_Hz    high-pass cutoff applied to detection signal
%                       (0 disables filtering; default: 0)
%
%     OPT.contrast
%       .Delta_bins     half-width of harmonic integration band [PSD bins]
%       .bg_bins        inner/outer bounds of background annulus [bins]
%
% -------------------------------------------------------------------------
% OUTPUTS
%
%   out   struct with fields:
%
%     out.meta
%       Fs             sampling frequency
%       T              total record duration
%       rpm            impeller speed
%       nblades        number of blades
%       fShaft         shaft frequency
%       fBPF           blade-passing frequency
%
%     out.PL
%       Upl, Vpl, Wpl  phase-locked velocity components
%       kept_lines     retained harmonic frequencies
%       leakage        per-harmonic energy capture diagnostics
%
%     out.residual
%       Ures, Vres, Wres   residual velocity components
%       Urms, Vrms, Wrms   RMS of residual components
%
%     out.triad
%       Emean          mean-flow kinetic energy
%       Epl            phase-locked kinetic energy
%       Eres           residual kinetic energy
%       frac           normalized energy fractions
%
%     out.figs_dir     directory where figures were saved (if enabled)
%
% -------------------------------------------------------------------------
% ASSUMPTIONS AND NOTES
%
%   • The velocity field is fully three-dimensional (U, V, W required).
%   • The time vector must be uniformly sampled.
%   • Harmonic selection is based solely on spectral contrast and does not
%     rely on phase averaging or blade-angle binning.
%   • Spectral filtering and contrast-band parameters affect harmonic
%     detection only; once harmonics are selected, the regression and energy
%     budgets are independent of these choices.
%
% -------------------------------------------------------------------------
% REFERENCES
%
%   If you use this function in published work, please cite the associated
%   manuscript and the Zenodo archive corresponding to this code release.
%
%% -------------------------------------------------------------------------

if ~exist('OPT','var'), OPT = struct; end
opt = OPT;

%% ------------ Defaults ------------
def.compute = struct( ...
    'triad_budget', true, ...
    'save_figs',    true, ...
    'show',         true);

def.save_dir = 'pod_figs';
def.dpi      = 200;

% harmonics:
%   Nmax       : max harmonic index (n) so that n*fShaft <= Nyq (default auto)
%   min_frac   : minimum contrast fraction to keep
%   guard_Hz   : half-band around fShaft for shaft removal when keep_shaft=false
%   keep_shaft : if true (default), shaft can be kept; if false, it is removed
def.harmonics = struct('Nmax',[],'min_frac',0.005,'guard_Hz',0.15,'keep_shaft',true);

% Welch PSD options (public knobs)
def.welch = struct( ...
    'Tw',   [], ...    % window length [s], default heuristic below
    'ovlp', 0.5, ...   % overlap fraction
    'nfft', [], ...    % empty => 2^nextpow2(Lw)
    'minL', 256);      % minimum segment length in samples

% PSD conditioning (public knobs)
def.psd = struct( ...
    'highpass_Hz',   0, ... % set 0 to disable
    'highpass_order', 4);

% Contrast bandwidths (public knobs, in PSD bins)
def.contrast = struct( ...
    'Delta_bins', 2, ...      % +/- bins around the line
    'bg_bins',    [3 8]);     % [W1 W2] bins for flanking annuli

opt           = setdefaults(opt, def);
opt.compute   = setdefaults(opt.compute, def.compute);
opt.harmonics = setdefaults(opt.harmonics, def.harmonics);
opt.welch     = setdefaults(getstruct_default(opt,'welch'), def.welch);
opt.psd       = setdefaults(getstruct_default(opt,'psd'), def.psd);
opt.contrast  = setdefaults(getstruct_default(opt,'contrast'), def.contrast);

% Backward-compatibility mapping
if isfield(opt,'welch_Tw') && isempty(opt.welch.Tw),   opt.welch.Tw = opt.welch_Tw;   end
if isfield(opt,'welch_ovlp') && isempty(opt.welch.ovlp), opt.welch.ovlp = opt.welch_ovlp; end

%% ------------ Unpack & normalize ------------
t    = DATA.t(:);
T    = numel(t);
Fs   = 1/mean(diff(t));
Ttot = t(end) - t(1);

rpm    = META.rpm;
nbl    = META.nblades;
fShaft = rpm/60;
fBPF   = nbl*fShaft;
Nyq    = Fs/2;

X = DATA.X;
Z = DATA.Z;
U = DATA.U;
W = DATA.W;

if ~isfield(DATA,'V') || isempty(DATA.V)
    error('DATA.V is required (3D velocity assumed).');
end
V = DATA.V;

if isfield(DATA,'S') && ~isempty(DATA.S)
    S = DATA.S;
else
    S = [];
end

% Normalize shapes to [Nz x Nx x T]; X,Z as [Nz x Nx]; S broadcastable
[U, W, V, X, Z, S] = normalize_shapes(U, W, V, X, Z, S);
[Nz, Nx, ~] = size(U);

% Masked time-means
Um = mean_masked(U,S);             % [Nz x Nx]
Vm = mean_masked(V,S);
Wm = mean_masked(W,S);

% Demeaned signals
U0 = U - Um;
V0 = V - Vm;
W0 = W - Wm;

% Axes (robust to orientation; match Nx, Nz)
xvec = linspace(min(X(:)), max(X(:)), Nx);
zvec = linspace(min(Z(:)), max(Z(:)), Nz).';

% Prepare folder only if we actually intend to save figures
if opt.compute.save_figs && ~isempty(opt.save_dir)
    if ~exist(opt.save_dir,'dir')
        mkdir(opt.save_dir);
    end
else
    opt.save_dir = '';  % ensure not used accidentally
end

% Welch parameters (default heuristic)
if isempty(opt.welch.Tw)
    % Prefer shaft-aware windowing: aim for ~20 shaft periods, capped.
    % Falls back safely even if fShaft is small.
    Tw0 = 20 / max(fShaft, 1e-6);
    opt.welch.Tw = min(6.0, max(1.5, Tw0));
end
Lw   = max(opt.welch.minL, round(opt.welch.Tw*Fs));
Nov  = round(opt.welch.ovlp*Lw);
if isempty(opt.welch.nfft)
    nfft = 2^nextpow2(Lw);
else
    nfft = opt.welch.nfft;
end

% ------------------------------------------------------------------
%% Phase-locked regression (shaft + BPF harmonics on shaft grid)
% ------------------------------------------------------------------

% Candidate lines: n*fShaft up to Nmax or Nyquist
if isempty(opt.harmonics.Nmax)
    Nmax_est = floor(min(Nyq / fShaft, nbl*10)); % up to ~10*BPF but ≤ Nyquist
else
    Nmax_est = opt.harmonics.Nmax;
end
fgrid = (2:Nmax_est)*fShaft;
fgrid = fgrid(fgrid < Nyq);

% PSD of plane-averaged kinetic proxy q(t)
[fen, Pen] = line_scan_q(U0, V0, W0, S, Fs, Lw, Nov, nfft, opt.psd);

% Background-subtracted contrast per candidate line (bin-based widths)
Cfrac = line_contrast(Pen, fen, fgrid, opt.contrast);

% Contrast-based selector
min_contrast = opt.harmonics.min_frac;
% legacy support if someone used a different field name
if isfield(opt.harmonics,'min_contrast_frac') && ~isempty(opt.harmonics.min_contrast_frac)
    min_contrast = opt.harmonics.min_contrast_frac;
end
kept = fgrid(Cfrac >= min_contrast).';

% Always keep BPF; keep shaft depending on keep_shaft flag
if opt.harmonics.keep_shaft
    must = unique([fShaft, fBPF]);
else
    must = fBPF;   % only BPF
end
must_keep = must(must < Nyq);
kept      = unique([kept(:); must_keep(:)]).';  % merge in column space

% Optional: force additional lines (Hz or integer multiples of fShaft)
if isfield(opt.harmonics,'force_lines') && ~isempty(opt.harmonics.force_lines)
    fl = opt.harmonics.force_lines(:).';
    if all(abs(fl - round(fl)) < 1e-6)     % integers => multiples of fShaft
        fl = fl * fShaft;
    end
    fl   = fl(fl < Nyq);
    kept = unique([kept(:); fl(:)]).';     % merge in column space
end

% If keep_shaft is false, remove anything near shaft and drop n=1
if ~opt.harmonics.keep_shaft
    kept = kept( abs(kept - fShaft) > opt.harmonics.guard_Hz );
    red  = round(kept / fShaft);
    kept = kept(red ~= 1);
end

fprintf('Kept %d harmonics. Reduced indices:\n', numel(kept));
if ~isempty(kept)
    disp(kept ./ fShaft);
end

% Design matrix for regression (cos/sin only; no DC)
Xreg = ones(T, 2*numel(kept));
for i = 1:numel(kept)
    Xreg(:, 2*i-1) = cos(2*pi*kept(i)*t);
    Xreg(:, 2*i  ) = sin(2*pi*kept(i)*t);
end

% Per-pixel masked LS regression
[Upl, Vpl, Wpl, leak] = fit_phase_locked(U0, V0, W0, S, Xreg, kept);

% Diagnostic figure: harmonic selection
if opt.compute.save_figs
    red_all = fgrid / fShaft;   % candidate grid
    yA      = Cfrac;            % normalized contrast (sum = 1)

    fh = figure('Color','w','Visible',onoff(opt.compute.show),'Name','phase_lines');
    grid on
    bar(red_all, yA, 0.9);
    xlabel('Reduced frequency (f/f_{shaft})');
    ylabel('Line energy fraction (contrast)');
    title(sprintf('%s — harmonic selection', META.name));
    export_fig_local(fh, opt, 'phase_lines');
end

% Residuals (all zero-mean by construction)
Ures = U0 - Upl;
Vres = V0 - Vpl;
Wres = W0 - Wpl;

% Basic variance diagnostics
Ur = rms3(Ures);
Vr = rms3(Vres);
Wr = rms3(Wres);
fprintf('RMS stats: Ures[max=%g min=%g], Vres[max=%g min=%g], Wres[max=%g min=%g]\n', ...
    max(Ur,[],'all'), min(Ur,[],'all'), ...
    max(Vr,[],'all'), min(Vr,[],'all'), ...
    max(Wr,[],'all'), min(Wr,[],'all'));

% NaN counts
nNanU = sum(~isfinite(Ures(:)));
nNanV = sum(~isfinite(Vres(:)));
nNanW = sum(~isfinite(Wres(:)));
fprintf('NaNs: Ures=%d, Vres=%d, Wres=%d\n', nNanU, nNanV, nNanW);

% Replace non-finite residuals with 0 (these regions had empty masks)
Ures(~isfinite(Ures)) = 0;
Vres(~isfinite(Vres)) = 0;
Wres(~isfinite(Wres)) = 0;

% ------------------------------------------------------------------
%% Energy budgets (mean / phase-locked / residual)
% ------------------------------------------------------------------
occ  = mean(S, 3, 'omitnan');             % 0..1 fraction of valid samples
num  = sum( (Um.^2 + Vm.^2 + Wm.^2) .* occ, 'all', 'omitnan');
den  = sum(occ, 'all', 'omitnan') + eps;
Emean = num / den;

EPL  = mean( (Upl.^2 + Vpl.^2 + Wpl.^2) .* S, 'all', 'omitnan');
ERES = mean( (Ures.^2 + Vres.^2 + Wres.^2) .* S, 'all', 'omitnan');

Ssum = Emean + EPL + ERES + eps;
frac = [Emean, EPL, ERES] / Ssum;

% RMS figures (PL & RES)
if opt.compute.triad_budget && opt.compute.save_figs
    % PL RMS
    fh = figure('Color','w','Visible',onoff(opt.compute.show),'Name','pl_rms');
    tiledlayout(1, 3,'Padding','compact','TileSpacing','compact');
    imagesc_axis(xvec,zvec, rms_masked(Upl,S), 'U_{PL} RMS');
    imagesc_axis(xvec,zvec, rms_masked(Vpl,S), 'V_{PL} RMS');
    imagesc_axis(xvec,zvec, rms_masked(Wpl,S), 'W_{PL} RMS');
    sgtitle(sprintf('%s — Phase-locked RMS', META.name));
    export_fig_local(fh, opt, 'pl_rms');

    % RES RMS
    fh = figure('Color','w','Visible',onoff(opt.compute.show),'Name','res_rms');
    tiledlayout(1, 3,'Padding','compact','TileSpacing','compact');
    imagesc_axis(xvec,zvec, rms_masked(Ures,S), 'U_{RES} RMS');
    imagesc_axis(xvec,zvec, rms_masked(Vres,S), 'V_{RES} RMS');
    imagesc_axis(xvec,zvec, rms_masked(Wres,S), 'W_{RES} RMS');
    sgtitle(sprintf('%s — Residual RMS', META.name));
    export_fig_local(fh, opt, 'res_rms');

    % Energy bars / pie
    fh = figure('Color','w','Visible',onoff(opt.compute.show),'Name','triad_energy');
    subplot(1,2,1);
    bar([Emean EPL ERES]);
    set(gca,'XTickLabel',{'Mean','PL','RES'});
    ylabel('Absolute energy (plane)'); grid on
    subplot(1,2,2);
    pie(frac, {'Mean','PL','RES'});
    sgtitle(sprintf('%s — Energy budgets', META.name));
    export_fig_local(fh, opt, 'triad_energy');
end

%% ----------------- Pack output -----------------
out = struct();
out.meta = struct('Fs',Fs,'T',Ttot,'rpm',rpm,'nblades',nbl, ...
    'fShaft',fShaft,'fBPF',fBPF, ...
    'welch_Tw',opt.welch.Tw,'welch_ovlp',opt.welch.ovlp, ...
    'psd_highpass_Hz',opt.psd.highpass_Hz, ...
    'contrast_Delta_bins',opt.contrast.Delta_bins, ...
    'contrast_bg_bins',opt.contrast.bg_bins);

out.PL = struct('Upl',Upl,'Vpl',Vpl,'Wpl',Wpl, ...
    'kept_lines',kept,'leakage',leak);

out.residual = struct('Ures',Ures,'Vres',Vres,'Wres',Wres, ...
    'Urms',rms3(Ures),'Vrms',rms3(Vres),'Wrms',rms3(Wres));

out.triad = struct('Emean',Emean,'Epl',EPL,'Eres',ERES,'frac',frac);

out.figs_dir = opt.save_dir;

if ~opt.compute.show
    close all;
end
end

%% ===================== Helpers =====================
function opt = setdefaults(opt, def)
fn = fieldnames(def);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(opt,f) || isempty(opt.(f))
        opt.(f) = def.(f);
    end
end
end

function s = getstruct_default(opt, f)
if isfield(opt,f) && isstruct(opt.(f))
    s = opt.(f);
else
    s = struct();
end
end

function [U, W, V, X, Z, S] = normalize_shapes(U, W, V, X, Z, S)
% Ensure [Nz x Nx x T] for fields, [Nz x Nx] for grids, mask broadcastable
[sz1, sz2, T] = size(U);

% If grids are transposed relative to fields, fix them
if ~isequal(size(X),[sz1 sz2]) && isequal(size(X),[sz2 sz1])
    X = X.';
    Z = Z.';
end

% If fields are transposed relative to grids, fix them
if ~isequal(size(U,1), size(X,1))
    U = permute(U,[2 1 3]);
    W = permute(W,[2 1 3]);
    if ~isempty(V), V = permute(V,[2 1 3]); end
end

% Default mask: all ones
if isempty(S), S = ones(size(U)); end
if ismatrix(S), S = repmat(S,1,1,T); end
end

function Y = mean_masked(F,S)
num = sum(F.*S,3,'omitnan');
den = sum(S,3,'omitnan');
den(den==0) = NaN;
Y = num./den;
end

function [fen, Pen] = line_scan_q(U,V,W,S,Fs,Lw,Nov,nfft, psdopt)
% PSD of plane-averaged kinetic proxy q(t)
q = 0.5 .* squeeze(mean(mean( (U.^2 + V.^2 + W.^2) .* S, ...
                              1,'omitnan'), 2,'omitnan'));
q = detrend(q,'linear');

% Optional high-pass to suppress ultra-low-frequency drift
f_hp = psdopt.highpass_Hz;
if ~isempty(f_hp) && f_hp > 0
    if f_hp >= (Fs/2)
        error('psd.highpass_Hz must be < Nyquist (Fs/2).');
    end
    ord = psdopt.highpass_order;
    [b,a] = butter(ord, f_hp/(Fs/2), 'high');
    q     = filtfilt(b,a,double(q));
end

[Pen, fen] = pwelch(q, hann(Lw,'periodic'), Nov, nfft, Fs, 'onesided');
end

function Cfrac = line_contrast(Pxx, f, fgrid, copt)
% Excess (background-subtracted) band power at each candidate line.
% Band widths defined in PSD bins, not in 1/Ttot.

if numel(f) < 2
    error('PSD frequency vector too short.');
end
df_psd = f(2) - f(1);

Delta = copt.Delta_bins * df_psd;
W1    = copt.bg_bins(1) * df_psd;
W2    = copt.bg_bins(2) * df_psd;

Ck = zeros(size(fgrid));
for i = 1:numel(fgrid)
    fk = fgrid(i);

    % Line band
    f1    = max(f(2), fk - Delta);
    f2    = min(f(end), fk + Delta);
    Pband = bandpower(Pxx, f, [f1 f2], 'psd');

    % Background from flanking annuli (median PSD * band width)
    L1   = [max(f(2), fk - W2), max(f(2), fk - W1)];
    L2   = [min(f(end), fk + W1), min(f(end), fk + W2)];
    mask = ((f>=L1(1) & f<=L1(2)) | (f>=L2(1) & f<=L2(2)));
    if any(mask)
        Pmed = median(Pxx(mask));
    else
        Pmed = median(Pxx); % fallback (rare; only near edges)
    end
    Ck(i) = max(0, Pband - Pmed * (f2 - f1));
end

Csum  = sum(Ck) + eps;
Cfrac = Ck / Csum;   % normalized contrast (for bars/thresholds)
end

function v = getfield_default(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function [Upl, Vpl, Wpl, leak] = fit_phase_locked(U,V,W,S,Xreg,kept)
[Nz,Nx,T] = size(U);
Upl = zeros(Nz,Nx,T);
Vpl = zeros(Nz,Nx,T);
Wpl = zeros(Nz,Nx,T);

% Coefficient-based leakage accumulators (combined over components)
npairs      = size(Xreg,2)/2;
Epl_accum   = zeros(npairs,1);   % mean-square per (cos,sin) pair
W_accum     = zeros(npairs,1);   % weights for normalization
Eres_accum  = 0;                 % residual mean-square (all comps)
Wres_accum  = 0;

% Per-pixel masked normal equations
for iz = 1:Nz
    for ix = 1:Nx
        s = squeeze(S(iz,ix,:)) > 0;
        if ~any(s), continue; end

        Xloc = Xreg(s,:);
        XtX  = Xloc.'*Xloc;
        lam  = 1e-10*trace(XtX)/size(XtX,1);
        G    = (XtX + lam*eye(size(XtX))) \ (Xloc.');

        occ         = mean(s);          % 0..1 time coverage for this pixel
        Eline_pixel = zeros(npairs,1);  % this pixel’s per-pair energy
        Eres_pix    = 0;                % this pixel’s residual energy

        % U component
        y = double(squeeze(U(iz,ix,s)));
        c = G*y;
        Upl(iz,ix,:) = Xreg*c;
        Eres_pix     = Eres_pix    + mean( (y - Xloc*c).^2 );
        Eline_pixel  = Eline_pixel + 0.5*(c(1:2:end).^2 + c(2:2:end).^2);

        % V component
        y = double(squeeze(V(iz,ix,s)));
        c = G*y;
        Vpl(iz,ix,:) = Xreg*c;
        Eres_pix     = Eres_pix    + mean( (y - Xloc*c).^2 );
        Eline_pixel  = Eline_pixel + 0.5*(c(1:2:end).^2 + c(2:2:end).^2);

        % W component
        y = double(squeeze(W(iz,ix,s)));
        c = G*y;
        Wpl(iz,ix,:) = Xreg*c;
        Eres_pix     = Eres_pix    + mean( (y - Xloc*c).^2 );
        Eline_pixel  = Eline_pixel + 0.5*(c(1:2:end).^2 + c(2:2:end).^2);

        % Accumulate plane-averaged, time-weighted energies
        Epl_accum  = Epl_accum  + occ * Eline_pixel;
        W_accum    = W_accum    + occ;
        Eres_accum = Eres_accum + occ * Eres_pix;
        Wres_accum = Wres_accum + occ;
    end
end

% Coefficient-based leakage metrics
Epl_line = Epl_accum ./ (W_accum + eps);
Eres     = Eres_accum ./ (Wres_accum + eps);
eta      = Epl_line ./ (Epl_line + Eres + eps);

% One entry per harmonic in `kept`
leak = struct('f', kept(:), ...
              'Epl_line',   Epl_line(:), ...
              'Eres_plane', Eres, ...
              'capture_eta',eta(:));
end

function r = rms3(F)
r = sqrt(mean(F.^2,3));
end

function v = onoff(tf)
if tf, v='on'; else, v='off'; end
end

function export_fig_local(fh, opt, tag)
% Save as <save_dir>/<tag>.[png|fig]
base = fullfile(opt.save_dir, tag);
print(fh, [base '.png'], '-dpng', sprintf('-r%d', getfield_default(opt,'dpi',200)));
visOld = get(fh,'Visible'); set(fh,'Visible','on');
try
    savefig(fh, [base '.fig']);
catch
    saveas(fh, [base '.fig']);
end
set(fh,'Visible',visOld);
end

function imagesc_axis(x,z,M,ttl)
% Strict: require axes to match matrix size (no silent repair)
Nz = size(M,1);
Nx = size(M,2);
if numel(x) ~= Nx || numel(z) ~= Nz
    error('imagesc_axis: axis sizes mismatch. Need numel(x)=%d, numel(z)=%d, got %d and %d.', ...
        Nx, Nz, numel(x), numel(z));
end

nexttile;
imagesc(x, z, M);
axis image; set(gca,'YDir','normal');
colorbar; xlabel('x'); ylabel('z'); title(ttl);
end

function R = rms_masked(F,S)
% RMS over time with mask S (assumed 0/1), returning [Nz x Nx]
num = sum((F.^2).*S,3,'omitnan');
den = sum(S,3,'omitnan');
den(den==0) = NaN;
R = sqrt(num./den);
end