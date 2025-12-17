function S = buildSMask(s, X, Z, freeSurf)
%BUILDSMASK Build a spatio-temporal validity mask for an x–z plane field.
%
%   S = buildSMask(s, X, Z, freeSurf) returns a mask S with the same size as s
%   (or broadcastable to it), where:
%       S = 1  -> valid sample
%       S = 0  -> excluded sample
%
%   Conventions (canonical layout):
%       X, Z : [Nz x Nx] coordinate grids
%       s    : [Nz x Nx x T] base mask/indicator (typically 0/1 or logical)
%       S    : [Nz x Nx x T] output mask (logical)
%
%   The mask is constructed by applying:
%       (1) a fixed polygon exclusion region in the (|x|, z) plane
%       (2) an optional free-surface exclusion polygon in the (x, z) plane
%
%   Inputs:
%       s        base indicator/mask, same grid as X,Z and time-resolved
%       X, Z     plane coordinate grids
%       freeSurf optional struct with fields:
%                   freeSurf.Points_0 : x coordinates of free-surface curve
%                   freeSurf.Points_2 : z coordinates of free-surface curve
%               Pass [] to disable free-surface masking.
%
%   Output:
%       S        logical mask [Nz x Nx x T]
%

%% ---- Parameter: fixed exclusion polygon in (|x|, z) ----
r_poly = [0.0617, 0.0586, 0.0263, 0.0230, 0.0201, 0.0000, 0.0000, 0.0454, 0.0626];
z_poly = [0.0594, 0.0627, 0.0627, 0.0595, 0.0200, 0.0200, 0.0040, 0.0040, 0.0220];

% Validate minimal sizes early (fail loudly if data are inconsistent)
if ndims(s) ~= 3
    error('buildSMask: s must be [Nz x Nx x T].');
end
[Nz, Nx, T] = size(s);
if ~isequal(size(X), [Nz Nx]) || ~isequal(size(Z), [Nz Nx])
    error('buildSMask: X and Z must be [Nz x Nx] and consistent with s.');
end

%% ---- Base mask ----
% Force to logical and preallocate output
S = logical(s);

%% ---- Exclusion 1: fixed polygon in (|x|, z) ----
in_fixed = inpolygon(abs(X), Z, r_poly, z_poly);

%% ---- Exclusion 2: optional free-surface polygon in (x, z) ----
use_free = ~isempty(freeSurf);
if use_free
    if ~isfield(freeSurf,'Points_0') || ~isfield(freeSurf,'Points_2')
        error('buildSMask: freeSurf must contain fields Points_0 and Points_2.');
    end

    xsurf = freeSurf.Points_0(:);
    zsurf = freeSurf.Points_2(:);

    % Sort by x to avoid self-intersections in polygon construction
    [xsurf, idx] = sort(xsurf);
    zsurf = zsurf(idx);

    % Close polygon by extending to a bounding box above the free surface
    x_ext = [xsurf;  0.25;  0.25; -0.25; -0.25];
    z_ext = [zsurf; zsurf(end); 1; 1; zsurf(1)];

    in_free = inpolygon(X, Z, x_ext, z_ext);
else
    in_free = false(Nz, Nx);
end

%% ---- Apply exclusions over time ----
% Exclusions are purely spatial; apply the same mask to each time slice.
excl = in_fixed | in_free;
for k = 1:T
    St = S(:,:,k);
    St(excl) = false;
    S(:,:,k) = St;
end
end
