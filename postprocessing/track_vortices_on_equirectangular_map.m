%% Post-process psi_a_vs_t: planisphere density + vortex trajectories
clear
close all
clc

%% User parameters
q_values = 1:10;

theta_cut = 0.18*pi;          % exclude polar caps
density_smoothing_sigma = 1.0;
min_distance_pixels = 8;

T_slope_ref = 0.05;           % barrier switch-off time

make_figures_visible = 'off';

arrow_length    = 0.10;       % arrow size in radians
arrow_linewidth = 3.8;

%% Paths
script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));

base_folder = fullfile(repository_root, 'output', 'Sweep_q');

output_root = fullfile(base_folder, 'Postprocess_vortex_tracks_planisphere');

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% Optional oslo colormap
try
    oslo_custom = Load_oslo_colormap();
    use_oslo = true;
catch
    use_oslo = false;
    warning(['Bundled Oslo colormap could not be loaded. ' ...
        'Falling back to the parula colormap.']);
end

%% Loop over q
for iq = 1:numel(q_values)

    q = q_values(iq);

    psi_folder = fullfile(base_folder, sprintf('q_%d', q), 'psi_a_vs_t');

    if ~isfolder(psi_folder)
        fprintf('Skipping q = %d: psi_a_vs_t not found\n', q)
        continue
    end

    files = dir(fullfile(psi_folder, 'psi_a_t_*.mat'));

    if isempty(files)
        fprintf('Skipping q = %d: no psi_a_t_*.mat files found\n', q)
        continue
    end

    [~, idx_sort] = sort({files.name});
    files = files(idx_sort);

    n_frames_all = numel(files);

    fprintf('Processing q = %d with %d saved snapshots\n', q, n_frames_all)

    q_output_folder = fullfile(output_root, sprintf('q_%d', q));

    if ~exist(q_output_folder, 'dir')
        mkdir(q_output_folder);
    end

    %% Load first frame for geometry
    S0 = load(fullfile(psi_folder, files(1).name), 'psi_a', 't', 'geom');

    geom = S0.geom;

    theta = geom.theta;
    phi   = geom.phi;

    if isvector(theta) && isvector(phi)
        [PHI, THETA] = meshgrid(phi, theta);
    else
        THETA = theta;
        PHI   = phi;
    end

    phi_vec   = PHI(1,:);
    theta_vec = THETA(:,1);

    polar_mask = (THETA < theta_cut) | (THETA > pi - theta_cut);

    %% Load all saved times and choose tracking start
    t_values_all = nan(n_frames_all,1);

    for it = 1:n_frames_all

        S_time = load(fullfile(psi_folder, files(it).name), 't');
        t_values_all(it) = S_time.t;

    end

    tracking_start_frame = find(t_values_all >= T_slope_ref, 1, 'first');

    if isempty(tracking_start_frame)

        warning(['q = %d: no frame found with t >= T_slope_ref. ' ...
            'Using first available frame.'], q)

        tracking_start_frame = 1;

    end

    tracking_files = files(tracking_start_frame:end);
    n_frames = numel(tracking_files);

    fprintf('  Tracking starts at frame %d, t = %.6f s\n', ...
        tracking_start_frame, t_values_all(tracking_start_frame))

    %% Number of vortices
    n_tracks = 2*q;

    vortex_theta = nan(n_frames, n_tracks);
    vortex_phi   = nan(n_frames, n_tracks);
    t_values     = nan(n_frames,1);

    previous_phi   = [];
    previous_theta = [];

    rho_plot = [];
    t_plot   = NaN;

    %% Detect and track vortices, starting after barrier switch-off
    for it = 1:n_frames

        S = load(fullfile(psi_folder, tracking_files(it).name), 'psi_a', 't');

        psi = S.psi_a;
        t_values(it) = S.t;

        rho = abs(psi).^2;

        if it == 1
            rho_plot = rho;
            t_plot = t_values(it);
        end

        rho_for_detection = rho;

        if density_smoothing_sigma > 0
            rho_for_detection = imgaussfilt(rho_for_detection, density_smoothing_sigma);
        end

        rho_for_detection(polar_mask) = max(rho_for_detection(:));

        candidates = find_density_holes_2d_periodic_phi( ...
            rho_for_detection, ...
            n_tracks, ...
            min_distance_pixels);

        cand_theta = THETA(sub2ind(size(THETA), candidates(:,1), candidates(:,2)));
        cand_phi   = PHI(sub2ind(size(PHI),   candidates(:,1), candidates(:,2)));

        cand_phi = mod(cand_phi, 2*pi);

        if it == 1

            [cand_phi, order] = sort(cand_phi);
            cand_theta = cand_theta(order);

            vortex_phi(it,:)   = cand_phi(:).';
            vortex_theta(it,:) = cand_theta(:).';

        else

            [matched_phi, matched_theta] = match_vortices_periodic_phi( ...
                previous_phi, previous_theta, ...
                cand_phi, cand_theta);

            vortex_phi(it,:)   = matched_phi(:).';
            vortex_theta(it,:) = matched_theta(:).';

        end

        previous_phi   = vortex_phi(it,:);
        previous_theta = vortex_theta(it,:);

    end

    %% Save tracking data
    output_mat = fullfile(q_output_folder, sprintf('vortex_tracks_q_%d.mat', q));

    save(output_mat, ...
        'q', ...
        't_values', ...
        't_values_all', ...
        'tracking_start_frame', ...
        'vortex_theta', ...
        'vortex_phi', ...
        'theta_cut', ...
        'T_slope_ref', ...
        'n_tracks', ...
        'rho_plot', ...
        't_plot', ...
        'geom');

    %% Plot planisphere
    fig = figure('Visible', make_figures_visible, ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',[3 3 24 13]);

    ax = axes(fig);
    ax.Position = [0.09 0.16 0.80 0.74];
    hold(ax,'on')

    imagesc(ax, phi_vec, theta_vec, rho_plot)

    % North pole at the top, South pole at the bottom.
    set(ax, 'YDir','reverse')

    if use_oslo
        colormap(ax, oslo_custom)
    else
        colormap(ax, parula)
    end

    cb = colorbar;
    cb.Ticks = [];
    cb.TickLabels = {};
    cb.Box = 'on';

    for n = 1:n_tracks

        %----------------------------------------------------------
        % Unwrapped version (used for arrow direction)
        %----------------------------------------------------------
        phi_unwrapped = unwrap(vortex_phi(:,n));

        %----------------------------------------------------------
        % Wrapped version (used for plotting)
        %----------------------------------------------------------
        phi_track   = mod(phi_unwrapped, 2*pi);
        theta_track = vortex_theta(:,n);

        %----------------------------------------------------------
        % Break trajectory when crossing phi = 0 / 2pi
        %----------------------------------------------------------
        jump_idx = find(abs(diff(phi_track)) > pi);

        phi_plot   = phi_track;
        theta_plot = theta_track;

        phi_plot(jump_idx+1)   = NaN;
        theta_plot(jump_idx+1) = NaN;

        %----------------------------------------------------------
        % Plot trajectory
        %----------------------------------------------------------
        hTrack = plot(ax, phi_plot, theta_plot, '-', ...
            'LineWidth', 3.8);

        track_color = hTrack.Color;

        %----------------------------------------------------------
        % Initial position
        %----------------------------------------------------------
        plot(ax, phi_track(1), theta_track(1), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k')

        %----------------------------------------------------------
        % Final arrow
        %----------------------------------------------------------
        [phi_end, theta_end, dphi_arrow, dtheta_arrow] = ...
            final_arrow_vector(phi_unwrapped, theta_track, arrow_length);

        phi_end = mod(phi_end, 2*pi);

        draw_custom_arrow_2d(ax, ...
            phi_end, theta_end, ...
            dphi_arrow, dtheta_arrow, ...
            track_color, ...
            arrow_linewidth, ...
            0.085, ...   % head length
            0.080)       % head width

    end

    xlabel(ax, '$\phi$', 'Interpreter','latex', 'FontSize', 24)
    ylabel(ax, '$\theta$', 'Interpreter','latex', 'FontSize', 24)

    title(ax, sprintf('Vortex tracking for $q=%d$', q), ...
        'Interpreter','latex', ...
        'FontSize', 24)

    xlim(ax, [0 2*pi])
    ylim(ax, [0 pi])

    xticks(ax, [0 pi/2 pi 3*pi/2 2*pi])
    xticklabels(ax, {'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})

    yticks(ax, [0 pi/4 pi/2 3*pi/4 pi])
    yticklabels(ax, {'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'})

    set(ax, ...
        'TickLabelInterpreter','latex', ...
        'FontSize', 22, ...
        'LineWidth', 1.5, ...
        'TickDir','out', ...
        'XColor','k', ...
        'YColor','k', ...
        'TickLength',[0.012 0.012], ...
        'Layer','top')

    box(ax,'on')

    drawnow

    %% Export
    exportgraphics(fig, ...
        fullfile(q_output_folder, sprintf('planisphere_tracks_q_%d.pdf', q)), ...
        'ContentType','vector')

    exportgraphics(fig, ...
        fullfile(q_output_folder, sprintf('planisphere_tracks_q_%d.png', q)), ...
        'Resolution',600)

    close(fig)

    fprintf('  Saved q = %d\n', q)

end

%% ========================================================================
% Local functions
%% ========================================================================

function candidates = find_density_holes_2d_periodic_phi(rho, n_holes, min_distance_pixels)

rho_work = rho;
candidates = nan(n_holes, 2);

n_theta = size(rho_work,1);
n_phi   = size(rho_work,2);

rho_max = max(rho_work(:));

for n = 1:n_holes

    [~, idx_min] = min(rho_work(:));
    [i_min, j_min] = ind2sub(size(rho_work), idx_min);

    candidates(n,:) = [i_min, j_min];

    i1 = max(1, i_min - min_distance_pixels);
    i2 = min(n_theta, i_min + min_distance_pixels);

    j_window = (j_min - min_distance_pixels):(j_min + min_distance_pixels);
    j_window = mod(j_window - 1, n_phi) + 1;

    rho_work(i1:i2, j_window) = rho_max;

end

end

function [matched_phi, matched_theta] = match_vortices_periodic_phi( ...
    previous_phi, previous_theta, cand_phi, cand_theta)

n_tracks = numel(previous_phi);

matched_phi   = nan(1,n_tracks);
matched_theta = nan(1,n_tracks);

available = true(numel(cand_phi),1);

for n = 1:n_tracks

    dphi = angular_distance_periodic(cand_phi, previous_phi(n));
    dtheta = cand_theta - previous_theta(n);

    dist2 = dtheta.^2 + dphi.^2;
    dist2(~available) = inf;

    [~, idx_best] = min(dist2);

    matched_phi(n)   = cand_phi(idx_best);
    matched_theta(n) = cand_theta(idx_best);

    available(idx_best) = false;

end

end

function dphi = angular_distance_periodic(phi1, phi2)

dphi = angle(exp(1i*(phi1 - phi2)));

end

function [phi_end, theta_end, dphi_arrow, dtheta_arrow] = ...
    final_arrow_vector(phi_track, theta_track, arrow_length)

valid = isfinite(phi_track) & isfinite(theta_track);

phi_track = phi_track(valid);
theta_track = theta_track(valid);

if numel(phi_track) < 2

    phi_end = phi_track(end);
    theta_end = theta_track(end);
    dphi_arrow = 0;
    dtheta_arrow = 0;
    return

end

phi_end = mod(phi_track(end), 2*pi);
theta_end = theta_track(end);

phi_prev = phi_track(end-1);
theta_prev = theta_track(end-1);

dphi = angular_distance_periodic(phi_track(end), phi_prev);
dtheta = theta_end - theta_prev;

norm_step = sqrt(dphi.^2 + dtheta.^2);

if norm_step < eps
    dphi_arrow = 0;
    dtheta_arrow = 0;
else
    dphi_arrow = arrow_length * dphi / norm_step;
    dtheta_arrow = arrow_length * dtheta / norm_step;
end

end

function draw_custom_arrow_2d(ax, x0, y0, dx, dy, color, lw, head_len, head_width)

    L = sqrt(dx^2 + dy^2);

    if L < eps
        return
    end

    ux = dx/L;
    uy = dy/L;

    % Arrow tip
    x_tip = x0 + dx;
    y_tip = y0 + dy;

    % End of shaft, before arrowhead
    x_base = x_tip - head_len*ux;
    y_base = y_tip - head_len*uy;

    % Perpendicular direction
    px = -uy;
    py = ux;

    % Shaft
    plot(ax, [x0 x_base], [y0 y_base], '-', ...
        'Color', color, ...
        'LineWidth', lw)

    % Arrowhead triangle
    x_head = [ ...
        x_tip, ...
        x_base + 0.5*head_width*px, ...
        x_base - 0.5*head_width*px];

    y_head = [ ...
        y_tip, ...
        y_base + 0.5*head_width*py, ...
        y_base - 0.5*head_width*py];

    patch(ax, x_head, y_head, color, ...
        'EdgeColor', color, ...
        'FaceColor', color)

end