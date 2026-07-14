%% Estimate sigma_* from highest-k Fourier mode of complex vortex displacement
clear
close all
clc

%% Parameters
q_values = 1:10;

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));

base_folder = fullfile(repository_root, 'output', 'Sweep_q');
track_root  = fullfile(base_folder, 'Postprocess_vortex_tracks_planisphere');

T_slope_ref = 0.05;

A_min_fit = 1e-9;
A_max_fit = 0.25;

min_fit_points = 8;
R2_min = 0.935;
segment_length_bonus = 0.002;

show_debug_logfits = true;

%% Output arrays
sigma_star             = nan(size(q_values));
sigma_star_err         = nan(size(q_values));
sigma_star_dimless     = nan(size(q_values));
sigma_star_dimless_err = nan(size(q_values));

sigma0          = nan(size(q_values));
sigma_ratio     = nan(size(q_values));
sigma_ratio_err = nan(size(q_values));

m_star         = nan(size(q_values));
A0_fit         = nan(size(q_values));
n_fit_pts      = nan(size(q_values));
R2_fit         = nan(size(q_values));
fit_start_time = nan(size(q_values));
fit_end_time   = nan(size(q_values));

%% Debug figure
if show_debug_logfits
    fig_log = figure('Color','w', ...
        'Units','centimeters', ...
        'Position',[2 2 35 18]);

    tl_log = tiledlayout(2,5, ...
        'TileSpacing','compact', ...
        'Padding','compact');
end

%% Loop over q
for iq = 1:numel(q_values)

    q = q_values(iq);

    track_file = fullfile(track_root, ...
        sprintf('q_%d', q), ...
        sprintf('vortex_tracks_q_%d.mat', q));

    if ~isfile(track_file)
        warning('Missing tracking file for q = %d', q)
        continue
    end

    S = load(track_file, ...
        't_values', ...
        'vortex_theta', ...
        'vortex_phi');

    t_all     = S.t_values(:);
    theta_all = S.vortex_theta;
    phi_all   = S.vortex_phi;

    %% Reference state: first state after barrier removal
    if min(t_all) < T_slope_ref && max(t_all) >= T_slope_ref
        ref_idx = find(t_all >= T_slope_ref, 1, 'first');
    else
        ref_idx = 1;
    end

    t     = t_all(ref_idx:end);
    theta = theta_all(ref_idx:end,:);
    phi   = phi_all(ref_idx:end,:);

    tau = t - t(1);

    n_frames = numel(t);
    n_tracks = size(phi,2);

    expected_tracks = 2*q;

    if n_tracks ~= expected_tracks
        warning('q = %d: expected %d vortices, found %d', ...
            q, expected_tracks, n_tracks)
    end

    %% Sort vortices according to post-barrier azimuthal order
    phi0_raw = mod(phi(1,:), 2*pi);

    [phi0, sort_idx] = sort(phi0_raw);

    phi   = phi(:, sort_idx);
    theta = theta(:, sort_idx);

    phi0   = phi0(:).';
    theta0 = theta(1,:);

    %% Unwrap phi trajectories consistently
    phi_unwrapped = nan(size(phi));

    for n = 1:n_tracks

        tmp = unwrap(phi(:,n));

        % Set the branch so that tmp(1) = phi0(n)
        tmp = tmp - tmp(1) + phi0(n);

        phi_unwrapped(:,n) = tmp;

    end

    %% Complex displacement from the post-barrier necklace
    % u_j(t) = delta_theta_j(t) + i delta_phi_j(t)

    delta_theta = theta - theta0;
    delta_phi   = phi_unwrapped - phi0;

    u = delta_theta + 1i*delta_phi;

    %% Highest-k mode supported by Nv = 2q vortices
    % m_max = Nv/2 = q

    m_theory = q;
    max_mode = floor(n_tracks/2);

    if m_theory > max_mode
        warning('q = %d: m_theory = %d exceeds max_mode = %d', ...
            q, m_theory, max_mode)
        continue
    end

    j = 0:(n_tracks-1);

    phase_factor = exp(-1i*2*pi*m_theory*j/n_tracks);

    c_m = nan(n_frames,1);

    for it = 1:n_frames
        c_m(it) = sum(u(it,:).*phase_factor)/n_tracks;
    end

    A = abs(c_m);
    logA = log(A);

    candidate = isfinite(A) & ...
        isfinite(logA) & ...
        A >= A_min_fit & ...
        A <= A_max_fit & ...
        tau > 0;

    if nnz(candidate) < min_fit_points
        warning('Not enough candidate points for q = %d, m = %d', ...
            q, m_theory)
        continue
    end

    %% Select best log-linear segment for fixed m = q
    [valid, p, R2_best, ~] = select_best_loglinear_segment( ...
        tau, ...
        logA, ...
        candidate, ...
        min_fit_points, ...
        R2_min, ...
        segment_length_bonus);

    if isempty(valid)
        warning('No sufficiently good fit found for q = %d, m = %d', ...
            q, m_theory)
        continue
    end

    x_fit = tau(valid);
    y_fit = logA(valid);

    m_star(iq) = m_theory;

    sigma_star(iq) = p(1);
    A0_fit(iq)     = exp(p(2));
    n_fit_pts(iq)  = numel(x_fit);
    R2_fit(iq)     = R2_best;

    fit_start_time(iq) = min(x_fit);
    fit_end_time(iq)   = max(x_fit);

    %% Fit uncertainty
    y_pred = polyval(p, x_fit);
    residuals = y_fit - y_pred;

    dof = n_fit_pts(iq) - 2;

    if dof > 0
        s2  = sum(residuals.^2) / dof;
        Sxx = sum((x_fit - mean(x_fit)).^2);
        sigma_star_err(iq) = sqrt(s2 / Sxx);
    else
        sigma_star_err(iq) = NaN;
    end

    %% Physical scales
    phys_file = fullfile(base_folder, ...
        sprintf('q_%d', q), ...
        'psi_a_vs_t', ...
        'psi_a_t_00001.mat');

    if ~isfile(phys_file)
        warning('Missing physical parameter file for q = %d', q)
        continue
    end

    S_phys = load(phys_file, 'hbar', 'm_a', 'geom');

    hbar = S_phys.hbar;
    m_a  = S_phys.m_a;
    R    = S_phys.geom.R_sphere;

    omega_unit = hbar/(m_a*R^2);

    kappa = 2*pi*hbar/m_a;
    a = 2*pi*R/(2*q);

    sigma0(iq) = pi*kappa/(4*a^2);

    sigma_star_dimless(iq)     = sigma_star(iq)/omega_unit;
    sigma_star_dimless_err(iq) = sigma_star_err(iq)/omega_unit;

    sigma_ratio(iq)     = sigma_star(iq)/sigma0(iq);
    sigma_ratio_err(iq) = sigma_star_err(iq)/sigma0(iq);

    fprintf(['q = %d: m = q = %d, sigma_* = %.4g +/- %.2g 1/s, ' ...
        'sigma_*/omega_unit = %.4g +/- %.2g, ' ...
        'sigma_*/sigma0 = %.4g +/- %.2g, ' ...
        'A0 = %.4g, points = %d, R2 = %.5f, ' ...
        'fit window = [%.4g, %.4g] s\n'], ...
        q, m_theory, ...
        sigma_star(iq), sigma_star_err(iq), ...
        sigma_star_dimless(iq), sigma_star_dimless_err(iq), ...
        sigma_ratio(iq), sigma_ratio_err(iq), ...
        A0_fit(iq), n_fit_pts(iq), R2_fit(iq), ...
        fit_start_time(iq), fit_end_time(iq));

    %% Debug log-fit plot
    if show_debug_logfits

        ax = nexttile(tl_log);
        hold(ax,'on')

        plot(ax, tau(candidate), logA(candidate), '.', ...
            'MarkerSize', 9)

        plot(ax, x_fit, y_fit, 'o', ...
            'MarkerSize', 4, ...
            'LineWidth', 1.0)

        tau_fit = linspace(min(x_fit), max(x_fit), 200);
        logA_fit = polyval(p, tau_fit);

        plot(ax, tau_fit, logA_fit, 'k-', ...
            'LineWidth', 2.0)

        xlabel(ax, '$t^\prime$ [s]', ...
            'Interpreter','latex')

        ylabel(ax, sprintf('$\\log |c_{%d}|$', m_theory), ...
            'Interpreter','latex')

        title(ax, sprintf('$q=%d$, $\\sigma_*=%.0f\\pm%.0f\\,\\mathrm{s}^{-1}$, $R^2=%.3f$', ...
            q, sigma_star(iq), sigma_star_err(iq), R2_fit(iq)), ...
            'Interpreter','latex', ...
            'FontSize',12)

        set(ax, ...
            'TickLabelInterpreter','latex', ...
            'FontSize',11, ...
            'LineWidth',1.0, ...
            'TickDir','out', ...
            'Layer','top', ...
            'XColor','k', ...
            'YColor','k')

        box(ax,'on')
        grid(ax,'on')

    end

end

%% Export debug log-fit figure

set(fig_log,'Renderer','painters')

exportgraphics(fig_log, ...
    'Sup_Mat_Log_fit.pdf', ...
    'ContentType','vector');

%% Results table
results_table = table(q_values(:), ...
    m_star(:), ...
    sigma_star(:), ...
    sigma_star_err(:), ...
    sigma_star_dimless(:), ...
    sigma_star_dimless_err(:), ...
    sigma0(:), ...
    sigma_ratio(:), ...
    sigma_ratio_err(:), ...
    A0_fit(:), ...
    n_fit_pts(:), ...
    R2_fit(:), ...
    fit_start_time(:), ...
    fit_end_time(:), ...
    'VariableNames', ...
    {'q', ...
    'm_star', ...
    'sigma_star', ...
    'sigma_star_err', ...
    'sigma_star_over_omega_unit', ...
    'sigma_star_over_omega_unit_err', ...
    'sigma0', ...
    'sigma_star_over_sigma0', ...
    'sigma_star_over_sigma0_err', ...
    'A0', ...
    'n_fit_points', ...
    'R2_fit', ...
    'fit_start_time', ...
    'fit_end_time'});

disp(results_table)

%% Final summary plot
fig_sigma = figure('Color','w', ...
    'Units','centimeters', ...
    'Position',[4 4 28 11]);

tl = tiledlayout(fig_sigma, 1, 2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% Optional point-vortex-model (PVM) reference curves
% These reference files are not part of the released repository. When
% present, they add a "Sphere PVM" comparison curve to the summary plots;
% otherwise that curve is silently omitted.
sphere_left_file  = fullfile(repository_root, 'sphere_linear_sigma_over_omega_unit.txt');
sphere_right_file = fullfile(repository_root, 'sphere_linear_sigma_over_sigma0.txt');

have_sphere_pvm = isfile(sphere_left_file) && isfile(sphere_right_file);

if have_sphere_pvm
    sphere_left  = readmatrix(sphere_left_file);
    sphere_right = readmatrix(sphere_right_file);
else
    warning(['Sphere PVM reference files not found:\n%s\n%s\n' ...
        'Skipping the "Sphere PVM" comparison curve.'], ...
        sphere_left_file, sphere_right_file);
end

%% Left panel

sigma0_dimless_theory = q_values.^2/2;


ax1 = nexttile(tl, 1);
hold(ax1,'on')
box(ax1,'on')
grid(ax1,'on')

errorbar(ax1, ...
    q_values, ...
    sigma_star_dimless, ...
    sigma_star_dimless_err, ...
    '-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k', ...
    'CapSize', 8)

plot(ax1, ...
    q_values, ...
    sigma0_dimless_theory, ...
    'k--', ...
    'LineWidth',2.5)

matlabRed = [0.8500 0.3250 0.0980];

if have_sphere_pvm
    plot(ax1, ...
        sphere_left(:,1), ...
        sphere_left(:,2), ...
        's-', ...
        'Color', matlabRed, ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor', matlabRed)
end

xlabel(ax1, '$q$', ...
    'Interpreter','latex', ...
    'FontSize', 26)

ylabel(ax1, '$\sigma_*/[\hbar/(m_aR^2)]$', ...
    'Interpreter','latex', ...
    'FontSize', 26)

set(ax1, ...
    'TickLabelInterpreter','latex', ...
    'FontSize', 21, ...
    'LineWidth', 1.4, ...
    'TickDir','out')

xlim(ax1, [min(q_values)-0.5, max(q_values)+0.5])
xticks(ax1, q_values)

legend_labels_1 = {'GPE', '$\sigma_0$'};
if have_sphere_pvm
    legend_labels_1{end+1} = 'Sphere PVM';
end

legend(ax1,...
    legend_labels_1,...
    'Interpreter','latex',...
    'Location','northwest')

%% Right panel
ax2 = nexttile(tl, 2);
hold(ax2,'on')
box(ax2,'on')
grid(ax2,'on')

errorbar(ax2, ...
    q_values, ...
    sigma_ratio, ...
    sigma_ratio_err, ...
    '-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k', ...
    'CapSize', 8)

xlabel(ax2, '$q$', ...
    'Interpreter','latex', ...
    'FontSize', 26)

ylabel(ax2, '$\sigma_*/\sigma_0$', ...
    'Interpreter','latex', ...
    'FontSize', 26)

set(ax2, ...
    'TickLabelInterpreter','latex', ...
    'FontSize', 21, ...
    'LineWidth', 1.4, ...
    'TickDir','out')

yline(ax2,1,'k--','LineWidth',2)

if have_sphere_pvm
    plot(ax2, ...
        sphere_right(:,1), ...
        sphere_right(:,2), ...
        's-', ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor','w')
end

ylim(ax2,[0 2])

xlim(ax2, [min(q_values)-0.5, max(q_values)+0.5])
xticks(ax2, q_values)

legend_labels_2 = {'GPE', '$\sigma_*/\sigma_0=1$'};
if have_sphere_pvm
    legend_labels_2{end+1} = 'Sphere PVM';
end

legend(ax2,...
    legend_labels_2,...
    'Interpreter','latex',...
    'Location','best')




%% ============================================================
% Robustness test of sigma*/sigma0
% ============================================================

Amax_list = [0.15 0.20 0.25 0.30];
R2_list   = [0.93 0.95 0.97];
Nmin_list = [6 8 10];

ratio_all = [];

for Amax_test = Amax_list

    for R2_test = R2_list

        for Nmin_test = Nmin_list

            sigma_ratio_tmp = nan(size(q_values));

            %--------------------------------------------------
            % Call the routine that extracts the growth rate sigma_*
            % using the same data but with
            %
            % A_max_fit     = Amax_test
            % R2_min        = R2_test
            % min_fit_points= Nmin_test
            %
            % and returns sigma_ratio_tmp(iq)
            %--------------------------------------------------

            sigma_ratio_tmp = compute_sigma_ratio_for_parameters( ...
                q_values, ...
                base_folder, ...
                track_root, ...
                T_slope_ref, ...
                A_min_fit, ...
                Amax_test, ...
                Nmin_test, ...
                R2_test, ...
                segment_length_bonus);

            ratio_all = [ratio_all ; sigma_ratio_tmp];

        end

    end

end

%% Statistics

ratio_mean = mean(ratio_all,1,'omitnan');
ratio_std  = std(ratio_all,0,1,'omitnan');

ratio_cv = 100*ratio_std./ratio_mean;

Trobust = table( ...
    q_values(:), ...
    ratio_mean(:), ...
    ratio_std(:), ...
    ratio_cv(:), ...
    'VariableNames', ...
    {'q','mean_ratio','std_ratio','percent_variation'});

disp(Trobust)

figure('Color','w')
hold on
box on
grid on

errorbar( ...
    q_values, ...
    ratio_mean, ...
    ratio_std, ...
    '-o', ...
    'LineWidth',2)

xlabel('$q$','Interpreter','latex')
ylabel('$\sigma_*/\sigma_0$','Interpreter','latex')

set(gca,...
    'TickLabelInterpreter','latex',...
    'FontSize',16)

%% =========================================================
% Export left panel only, with inset sigma*/sigma0
% =========================================================
matlabBlue = [0 0.4470 0.7410];
matlabRed  = [0.8500 0.3250 0.0980];

fig_left = figure( ...
    'Color','w', ...
    'Position',[100 100 700 550]);

ax = axes(fig_left);
hold(ax,'on')
box(ax,'on')
grid(ax,'on')

errorbar(ax,...
    q_values,...
    sigma_star_dimless,...
    sigma_star_dimless_err,...
    '-o',...
    'Color',matlabBlue,...
    'LineWidth',2.5,...
    'MarkerSize',8,...
    'MarkerFaceColor','w',...
    'MarkerEdgeColor',matlabBlue,...
    'CapSize',8)

plot(ax,...
    q_values,...
    sigma0_dimless_theory,...
    'k--',...
    'LineWidth',2.5)

matlabRed = [0.8500 0.3250 0.0980];

if have_sphere_pvm
    plot(ax, ...
        sphere_left(:,1), ...
        sphere_left(:,2), ...
        's-', ...
        'Color', matlabRed, ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor', matlabRed)
end

xlabel(ax,'$q$',...
    'Interpreter','latex',...
    'FontSize',26)

ylabel(ax,'$\sigma_*/[\hbar/(m_aR^2)]$',...
    'Interpreter','latex',...
    'FontSize',26)

set(ax,...
    'TickLabelInterpreter','latex',...
    'FontSize',21,...
    'LineWidth',1.4,...
    'TickDir','out')

xlim(ax,[min(q_values)-0.5,max(q_values)+0.5])
xticks(ax,q_values)

legend_labels_left = {'GPE', '$\sigma_0$'};
if have_sphere_pvm
    legend_labels_left{end+1} = 'Sphere PVM';
end

legend(ax,...
    legend_labels_left,...
    'Interpreter','latex',...
    'Location','southeast')

%% Inset: sigma*/sigma0 vs q

ax_in = axes(fig_left, ...
    'Position',[0.24 0.60 0.32 0.30]);   % top-left

hold(ax_in,'on')
box(ax_in,'on')
grid(ax_in,'on')

matlabBlue = [0 0.4470 0.7410];

errorbar(ax_in,...
    q_values,...
    sigma_ratio,...
    sigma_ratio_err,...
    'o-',...
    'Color',matlabBlue,...
    'LineWidth',1.8,...
    'MarkerSize',5.5,...
    'MarkerFaceColor','w',...
    'MarkerEdgeColor',matlabBlue,...
    'CapSize',5)

yline(ax_in,1,'k--','LineWidth',1.2)

if have_sphere_pvm
    plot(ax_in,...
        sphere_right(:,1),...
        sphere_right(:,2),...
        's-',...
        'Color',matlabRed,...
        'LineWidth',1.8,...
        'MarkerSize',5.5,...
        'MarkerFaceColor','w',...
        'MarkerEdgeColor',matlabRed)
end

xlabel(ax_in,'$q$',...
    'Interpreter','latex',...
    'FontSize',18)

ylabel(ax_in,'$\sigma_*/\sigma_0$',...
    'Interpreter','latex',...
    'FontSize',18)

set(ax_in,...
    'TickLabelInterpreter','latex',...
    'FontSize',16,...
    'LineWidth',1.2,...
    'TickDir','out',...
    'Layer','top')

xlim(ax_in,[min(q_values)-0.5,max(q_values)+0.5])
xticks(ax_in,2:2:max(q_values))

ylim(ax_in,[0 2])
yticks(ax_in,0:0.5:2)

%% Export

set(fig_left,'Renderer','painters')

exportgraphics(fig_left,...
    'sigma_star_vs_q_with_inset.pdf',...
    'ContentType','vector');

exportgraphics(fig_left,...
    'sigma_star_vs_q_with_inset.png',...
    'Resolution',600);

%% ========================================================================
% Local functions
%% ========================================================================

function [valid, best_p, best_R2, best_score] = select_best_loglinear_segment( ...
    tau, logA, candidate, min_fit_points, R2_min, segment_length_bonus)

valid = [];
best_p = [NaN NaN];
best_R2 = NaN;
best_score = -Inf;

idx_all = find(candidate(:));

if numel(idx_all) < min_fit_points
    return
end

breaks = find(diff(idx_all) > 1);

block_start = [1; breaks + 1];
block_end   = [breaks; numel(idx_all)];

for b = 1:numel(block_start)

    block_idx = idx_all(block_start(b):block_end(b));

    if numel(block_idx) < min_fit_points
        continue
    end

    for i1 = 1:(numel(block_idx)-min_fit_points+1)

        for i2 = (i1+min_fit_points-1):numel(block_idx)

            seg_idx = block_idx(i1:i2);

            x = tau(seg_idx);
            y = logA(seg_idx);

            if any(~isfinite(x)) || any(~isfinite(y))
                continue
            end

            p = polyfit(x, y, 1);

            if p(1) <= 0
                continue
            end

            y_fit = polyval(p, x);

            SS_res = sum((y - y_fit).^2);
            SS_tot = sum((y - mean(y)).^2);

            if SS_tot <= 0
                continue
            end

            R2 = 1 - SS_res/SS_tot;

            if R2 < R2_min
                continue
            end

            n_seg = numel(seg_idx);

            early_penalty = 0.08 * (min(x) / max(tau(candidate)));

            score = R2 + segment_length_bonus*n_seg - early_penalty;

            if score > best_score

                best_score = score;

                valid = false(size(candidate));
                valid(seg_idx) = true;

                best_p = p;
                best_R2 = R2;

            end

        end

    end

end

end

function sigma_ratio_tmp = compute_sigma_ratio_for_parameters( ...
    q_values, ...
    base_folder, ...
    track_root, ...
    T_slope_ref, ...
    A_min_fit, ...
    A_max_fit, ...
    min_fit_points, ...
    R2_min, ...
    segment_length_bonus)

sigma_ratio_tmp = nan(size(q_values));

for iq = 1:numel(q_values)

    q = q_values(iq);

    track_file = fullfile(track_root, ...
        sprintf('q_%d', q), ...
        sprintf('vortex_tracks_q_%d.mat', q));

    if ~isfile(track_file)
        continue
    end

    S = load(track_file, ...
        't_values', ...
        'vortex_theta', ...
        'vortex_phi');

    t_all     = S.t_values(:);
    theta_all = S.vortex_theta;
    phi_all   = S.vortex_phi;

    if min(t_all) < T_slope_ref && max(t_all) >= T_slope_ref
        ref_idx = find(t_all >= T_slope_ref, 1, 'first');
    else
        ref_idx = 1;
    end

    t     = t_all(ref_idx:end);
    theta = theta_all(ref_idx:end,:);
    phi   = phi_all(ref_idx:end,:);

    tau = t - t(1);

    n_tracks = size(phi,2);

    phi0_raw = mod(phi(1,:), 2*pi);

    [phi0, sort_idx] = sort(phi0_raw);

    phi   = phi(:, sort_idx);
    theta = theta(:, sort_idx);

    phi0   = phi0(:).';
    theta0 = theta(1,:);

    phi_unwrapped = nan(size(phi));

    for n = 1:n_tracks

        tmp = unwrap(phi(:,n));
        tmp = tmp - tmp(1) + phi0(n);

        phi_unwrapped(:,n) = tmp;

    end

    delta_theta = theta - theta0;
    delta_phi   = phi_unwrapped - phi0;

    u = delta_theta + 1i*delta_phi;

    m_theory = q;
    max_mode = floor(n_tracks/2);

    if m_theory > max_mode
        continue
    end

    j = 0:(n_tracks-1);

    phase_factor = exp(-1i*2*pi*m_theory*j/n_tracks);

    c_m = nan(numel(t),1);

    for it = 1:numel(t)
        c_m(it) = sum(u(it,:).*phase_factor)/n_tracks;
    end

    A = abs(c_m);
    logA = log(A);

    candidate = isfinite(A) & ...
        isfinite(logA) & ...
        A >= A_min_fit & ...
        A <= A_max_fit & ...
        tau > 0;

    if nnz(candidate) < min_fit_points
        continue
    end

    [valid, p, ~, ~] = select_best_loglinear_segment( ...
        tau, ...
        logA, ...
        candidate, ...
        min_fit_points, ...
        R2_min, ...
        segment_length_bonus);

    if isempty(valid)
        continue
    end

    sigma_star_tmp = p(1);

    phys_file = fullfile(base_folder, ...
        sprintf('q_%d', q), ...
        'psi_a_vs_t', ...
        'psi_a_t_00001.mat');

    if ~isfile(phys_file)
        continue
    end

    S_phys = load(phys_file, 'hbar', 'm_a', 'geom');

    hbar = S_phys.hbar;
    m_a  = S_phys.m_a;
    R    = S_phys.geom.R_sphere;

    kappa = 2*pi*hbar/m_a;
    a = 2*pi*R/(2*q);

    sigma0_tmp = pi*kappa/(4*a^2);

    sigma_ratio_tmp(iq) = sigma_star_tmp/sigma0_tmp;

end

end