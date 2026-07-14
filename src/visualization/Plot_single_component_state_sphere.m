function fig = Plot_single_component_state_sphere(psi, L_band, method, varargin)
%PLOT_FRAME_SINGLE_COMPONENT_SPHERE
% Plot density and phase of a single wavefunction on the sphere.
%
% INPUTS
%   psi     : wavefunction
%   L_band  : spherical harmonic band-limit
%   method  : SSHT sampling method (e.g. 'MW')
%
% OPTIONAL NAME-VALUE PAIRS
%   'Visible'     : 'on' (default) or 'off'
%   'FigureTitle' : title for the whole figure (default: '')
%
% OUTPUT
%   fig     : figure handle

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FigureTitle', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    fig_visible = char(p.Results.Visible);
    fig_title   = char(p.Results.FigureTitle);

    old_vis = get(groot, 'DefaultFigureVisible');
    cleanupObj = onCleanup(@() set(groot, 'DefaultFigureVisible', old_vis)); %#ok<NASGU>

    if strcmpi(fig_visible, 'off')
        set(groot, 'DefaultFigureVisible', 'off');
    end

    fig = figure('Visible', fig_visible, 'Position', [100, 100, 1200, 500]);

    tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    rho   = abs(psi).^2;
    phase = angle(psi);

    % -----------------------------
    % Left: density
    % -----------------------------
    ax1 = nexttile(tl, 1);
    axes(ax1);
    set(fig,'CurrentAxes',ax1)

    ssht_plot_sphere(rho, L_band, ...
        'Method', method, ...
        'Type', 'colour', ...
        'ColourBar', false);

    title('$|\psi|^2$', 'Interpreter', 'latex', 'FontSize', 18);
    colormap(ax1, gray);
    set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
    xlabel('$x/R$','Interpreter','latex','FontSize',15)
    ylabel('$y/R$','Interpreter','latex','FontSize',15)
    zlabel('$z/R$','Interpreter','latex','FontSize',15)

    % -----------------------------
    % Right: phase
    % -----------------------------
    ax2 = nexttile(tl, 2);
    axes(ax2);
    set(fig,'CurrentAxes',ax2)

    ssht_plot_sphere(phase, L_band, ...
        'Method', method, ...
        'Type', 'colour', ...
        'ColourBar', false);

    title('$\arg(\psi)$', 'Interpreter', 'latex', 'FontSize', 18);
    colormap(ax2, hsv);
    set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
    clim(ax2, [-pi pi]);
    xlabel('$x/R$','Interpreter','latex','FontSize',15)
    ylabel('$y/R$','Interpreter','latex','FontSize',15)
    zlabel('$z/R$','Interpreter','latex','FontSize',15)

    if ~isempty(fig_title)
        sgtitle(fig_title, 'Interpreter', 'latex', 'FontSize', 18);
    end

    set(fig, 'Visible', fig_visible);
    drawnow
end