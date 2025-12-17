%% run_demo.m — Example driver for flow_decomposition
%
% This script demonstrates the use of FLOW_DECOMPOSITION on a time-resolved
% velocity field extracted on an x–z plane.
%
% -------------------------------------------------------------------------
% DATA LAYOUT AND CONVENTIONS
%
% All spatial fields must follow the same canonical layout:
%
%   - Rows    correspond to the z-direction   (Nz points)
%   - Columns correspond to the x-direction   (Nx points)
%
% Required sizes:
%   X, Z        : [Nz x Nx]      coordinate grids
%   U, V, W     : [Nz x Nx x T]  velocity components
%   S (mask)    : [Nz x Nx] or [Nz x Nx x T]
%   t           : [T x 1]        uniformly sampled time vector
%
% With this convention:
%   - x-axis vector is obtained as X(1,:)
%   - z-axis vector is obtained as Z(:,1)
%
% The demo dataset provided with the Zenodo archive has already been
% canonicalized to satisfy these conventions.
%
% -------------------------------------------------------------------------

clear; close all; clc;

% Add source folder (adjust if repository layout differs)
addpath(fullfile('..','src'));

%% Load demo dataset
S0 = load(fullfile('..','data','dataset_demo.mat'));
ds = S0.dataset_demo;

%% Build spatial mask
% The mask must be consistent with the [Nz x Nx x T] data layout.
S = buildSMask(ds.s, ds.X, ds.Z, []);

%% Assemble DATA structure
DATA = struct();
DATA.t = ds.t(:);        % [T x 1]
DATA.X = ds.X;           % [Nz x Nx]
DATA.Z = ds.Z;           % [Nz x Nx]
DATA.U = ds.U;           % [Nz x Nx x T]
DATA.V = ds.V;           % [Nz x Nx x T]
DATA.W = ds.W;           % [Nz x Nx x T]
DATA.S = S;                   % [Nz x Nx] or [Nz x Nx x T]

%% Metadata
META = struct();
META.rpm     = 110;
META.nblades = 4;
META.name    = 'DES 19cm';

%% Options
OPT = struct();

% Output and diagnostics
OPT.save_dir = fullfile('..','Figures');
OPT.compute  = struct( ...
    'triad_budget', true, ...
    'save_figs',    true, ...
    'show',         true );

% Harmonic detection and selection
Ttot = DATA.t(end) - DATA.t(1);
OPT.harmonics = struct( ...
    'min_frac',    0.2, ...
    'keep_shaft',  false, ...
    'guard_Hz',    max(0.2, 5/max(Ttot,eps)), ...
    'force_lines', [] );

% Spectral analysis controls
OPT.welch = struct('Tw', min(5.0, 0.75*Ttot), 'ovlp', 0.5);
OPT.psd   = struct('highpass_Hz', 0);   % disabled by default

%% Run decomposition
out = flow_decomposition(DATA, META, OPT);

%% Post-processing: velocity magnitudes
U = DATA.U; V = DATA.V; W = DATA.W; S = DATA.S;

% Masked mean velocity
Um = masked_mean(U,S);
Vm = masked_mean(V,S);
Wm = masked_mean(W,S);

mag_mean = sqrt(Um.^2 + Vm.^2 + Wm.^2);

% RMS of phase-locked and residual components
mag_pl  = sqrt( ...
    masked_rms(out.PL.Upl,S).^2 + ...
    masked_rms(out.PL.Vpl,S).^2 + ...
    masked_rms(out.PL.Wpl,S).^2 );

mag_res = sqrt( ...
    masked_rms(out.residual.Ures,S).^2 + ...
    masked_rms(out.residual.Vres,S).^2 + ...
    masked_rms(out.residual.Wres,S).^2 );

%% Visualization
colorlim = [0 0.3];

plotModes(DATA.X(1,:), DATA.Z(:,1), mag_mean, colorlim);
title('Mean flow |U|');

plotModes(DATA.X(1,:), DATA.Z(:,1), mag_pl, colorlim);
title('Phase-locked RMS |U_{PL}|');

plotModes(DATA.X(1,:), DATA.Z(:,1), mag_res, colorlim);
title('Residual RMS |U_{RES}|');

%% ---------------- Local utilities ----------------
function M = masked_mean(F,S)
num = sum(F.*S, 3, 'omitnan');
den = sum(S,   3, 'omitnan');
den(den==0) = NaN;
M = num./den;
end

function R = masked_rms(F,S)
num = sum((F.^2).*S, 3, 'omitnan');
den = sum(S,        3, 'omitnan');
den(den==0) = NaN;
R = sqrt(num./den);
end
