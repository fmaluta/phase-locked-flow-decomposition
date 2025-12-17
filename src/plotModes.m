function h = plotModes(x, z, M, colorlim, opts)
%PLOTMODES Display a scalar field on an x–z plane.
%
%   h = plotModes(x, z, M, colorlim) displays matrix M on axes (x,z) using
%   imagesc and applies optional color limits.
%
%   Conventions:
%       x : [1 x Nx] or [Nx x 1]   x-coordinates (columns of M)
%       z : [Nz x 1] or [1 x Nz]   z-coordinates (rows of M)
%       M : [Nz x Nx]             scalar field
%
%   colorlim:
%       []        -> do not set clim
%       [cmin cmax] with cmin<cmax -> set clim([cmin cmax])
%
%   opts (optional struct):
%       .title        string   (default: '')
%       .xlabel       string   (default: 'x [m]')
%       .ylabel       string   (default: 'z [m]')
%       .cbar_label   string   (default: '')
%       .colormap     string   (default: '')  % '' -> keep current colormap
%       .overlay_poly logical  (default: true)
%       .poly_x       vector   % polygon x-vertices (default: built-in)
%       .poly_z       vector   % polygon z-vertices (default: built-in)
%       .poly_color   1x3      (default: [0.7 0.7 0.7])
%       .poly_lw      scalar   (default: 2)
%
%   Output:
%       h.fig, h.ax, h.im, h.cb, h.patch  graphics handles

if nargin < 4, colorlim = []; end
if nargin < 5, opts = struct(); end

% Defaults
def = struct();
def.title        = '';
def.xlabel       = 'x [m]';
def.ylabel       = 'z [m]';
def.cbar_label   = '';
def.colormap     = '';
def.overlay_poly = true;
def.poly_color   = [0.7 0.7 0.7];
def.poly_lw      = 2;

% Default polygon (same region used in buildSMask, but mirrored)
rpolyg_p = [0.0454,0.0626,0.0617,0.0586,0.0263,0.0230,0.0201];
zpolyg_p = [0.0040,0.0220,0.0594,0.0627,0.0627,0.0595,0.0200];
def.poly_x = [rpolyg_p, -fliplr(rpolyg_p)];
def.poly_z = [zpolyg_p,  fliplr(zpolyg_p)];

opts = setdefaults(opts, def);

% Force axis vectors as columns/rows for size checks
x = x(:).';     % 1 x Nx
z = z(:);       % Nz x 1

[Nz, Nx] = size(M);
if numel(x) ~= Nx || numel(z) ~= Nz
    error('plotModes: size mismatch. Need numel(x)=%d, numel(z)=%d for M[%d x %d].', ...
        Nx, Nz, Nz, Nx);
end

% Figure and axes
h.fig = figure('Color','w');
h.ax  = axes('Parent', h.fig);

h.im = imagesc(h.ax, x, z, M);
set(h.ax, 'YDir','normal');
axis(h.ax, 'image');
grid(h.ax, 'off');

% Color limits (only if valid)
if ~isempty(colorlim) && numel(colorlim)==2 && colorlim(2) > colorlim(1)
    clim(h.ax, colorlim);
end

% Colormap (only if requested)
if ~isempty(opts.colormap)
    colormap(h.ax, opts.colormap);
end

% Overlay polygon (optional)
h.patch = gobjects(1);
if opts.overlay_poly
    hold(h.ax, 'on');
    h.patch = fill(h.ax, opts.poly_x, opts.poly_z, opts.poly_color, ...
        'LineWidth', opts.poly_lw);
    hold(h.ax, 'off');
end

% Labels/title
xlabel(h.ax, opts.xlabel);
ylabel(h.ax, opts.ylabel);
if ~isempty(opts.title)
    title(h.ax, opts.title);
end

% Colorbar
h.cb = colorbar(h.ax);
if ~isempty(opts.cbar_label)
    h.cb.Label.String = opts.cbar_label;
end
end

%% ---- local helper ----
function opt = setdefaults(opt, def)
fn = fieldnames(def);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(opt,f) || isempty(opt.(f))
        opt.(f) = def.(f);
    end
end
end
