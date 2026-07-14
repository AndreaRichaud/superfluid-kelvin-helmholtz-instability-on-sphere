function Prepare_single_component_initial_state_sphere(q, params, output_folder)
%PREPARE_SINGLE_COMPONENT_INITIAL_STATE_SPHERE Prepare the initial state.
%
% Construct a single-component condensate on a spherical shell with:
%   - approximately uniform density in the two hemispheres;
%   - a high, static Gaussian barrier at the equator;
%   - two repulsive pinning potentials at the poles;
%   - a patchwise phase +q*phi in the northern hemisphere and -q*phi in
%     the southern hemisphere.
%
% This function does not lower the equatorial barrier in time. It prepares
% the state subsequently used by Run_single_component_dynamics_sphere.
%
% INPUTS
%   q             - integer vortex charge at the two poles
%   params        - structure containing numerical and physical parameters
%   output_folder - case output directory
%
% The function saves the prepared state, log file, intermediate frames,
% final-state figures, and imaginary-time diagnostics in output_folder.

    %-----------------------------%
    % Directory and logging setup

    if nargin < 2 || isempty(params)
        params = struct();
    end

    params = set_default_params(params);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    temp_img_folder = fullfile(output_folder, 'Temp_img');
    if ~exist(temp_img_folder, 'dir')
        mkdir(temp_img_folder);
    end

    log_file = fullfile(output_folder, 'log_compute_ground_state.txt');
    if exist(log_file, 'file')
        delete(log_file);
    end

    diary(log_file);
    diary on;

    t_start = tic;

    try
        fprintf('============================================================\n');
        fprintf('Prepare_single_component_initial_state_sphere started\n');
        fprintf('q = %d\n', q);
        fprintf('Output folder: %s\n', output_folder);
        fprintf('Start time: %s\n', datestr(now));
        fprintf('============================================================\n');

        %--------------------------------------------------------------%
        % Path setup
        %--------------------------------------------------------------%
        this_file_folder = fileparts(mfilename('fullpath'));
        addpath(genpath(this_file_folder));

        %--------------------------------------------------------------%
        % Physical constants
        %--------------------------------------------------------------%
        hbar = 1.0545718e-34;      % [J s]
        amu  = 1.66054e-27;        % [kg]
        a0   = 5.291777e-11;       % [m]

        %--------------------------------------------------------------%
        % Physical parameters: single component only
        %--------------------------------------------------------------%
        m_a = 23 * amu;            % Sodium-23
        N_a = params.N_a;

        R_sphere = params.R_sphere;   % [m]
        L_z      = params.L_z;        % [m]

        a_aa = params.a_aa * a0;      % [m]
        g_a  = 4*pi*hbar^2*a_aa/m_a;  % 3D coupling

        % effective quasi-2D coupling used only to estimate mu
        g2D_est = g_a / L_z;
        n2D_est = N_a / (4*pi*R_sphere^2);
        mu_est  = g2D_est * n2D_est;

        %--------------------------------------------------------------%
        % Numerical parameters
        %--------------------------------------------------------------%
        n_iterations     = params.n_iterations_imag;
        dt               = params.dt_imag;
        L_band           = params.L_band;
        sample_frequency = params.sample_frequency;

        %--------------------------------------------------------------%
        % Geometry
        %--------------------------------------------------------------%
        geom = Build_geometry_sphere(R_sphere, L_band, 'MW');

        fprintf('Geometry built successfully.\n');
        fprintf('L_band = %d\n', L_band);
        fprintf('n_iterations_imag = %d\n', n_iterations);
        fprintf('dt_imag = %.3e\n', dt);

        %--------------------------------------------------------------%
        % Static preparation potential:
        % equatorial barrier + optional pinning at poles
        %--------------------------------------------------------------%
        sigma_barrier = params.sigma_barrier;
        V0_barrier    = params.V0_over_mu * mu_est;

        if params.use_pinning_poles
            V0_pin   = params.V0_pin_over_mu * mu_est;
            sigma_pin = params.sigma_pin;
        else
            V0_pin   = 0;
            sigma_pin = params.sigma_pin;
        end

        Mat_V = build_preparation_potential_sphere( ...
            geom, V0_barrier, sigma_barrier, V0_pin, sigma_pin);

        fprintf('Preparation potential built.\n');
        fprintf('Estimated mu = %.3e J\n', mu_est);
        fprintf('Equatorial barrier height V0_barrier = %.3e J\n', V0_barrier);
        fprintf('Equatorial barrier width sigma_barrier = %.4f rad\n', sigma_barrier);
        fprintf('Polar pinning height V0_pin = %.3e J\n', V0_pin);
        fprintf('Polar pinning width sigma_pin = %.4f rad\n', sigma_pin);

        %--------------------------------------------------------------%
        % Initial condition: density dip at equator + patchwise phase
        %--------------------------------------------------------------%
        psi = build_initial_patchwise_state_sphere(geom, q, params);

        % normalize with geometry weights
        psi = psi / sqrt(sum(abs(psi(:)).^2 .* geom.W(:)));

        fprintf('Initial single-component state generated.\n');

        %--------------------------------------------------------------%
        % Diagnostics allocation
        %--------------------------------------------------------------%
        n_samples = ceil(n_iterations / sample_frequency) + 1;

        vec_mu  = zeros(n_samples,1);
        vec_Ene = zeros(n_samples,1);
        vec_t   = zeros(n_samples,1);

        i_sampling = 1;

        %--------------------------------------------------------------%
        % Imaginary-time evolution with STATIC preparation potential
        %--------------------------------------------------------------%
        for i = 1:n_iterations

            [psi, Hpsi] = Imaginary_time_evolve_sphere_single_component( ...
                psi, m_a, g_a, N_a, L_z, Mat_V, geom, hbar, dt);

            if i == 1 || mod(i, sample_frequency) == 0

                elapsed_time = toc(t_start);
                fprintf('Progress = %.2f %% in %.0f s\n', ...
                    100*i/n_iterations, elapsed_time);

                vec_mu(i_sampling) = real(sum(conj(psi(:)) .* Hpsi(:) .* geom.W(:)));

                vec_Ene(i_sampling) = Compute_total_energy_single_component_sphere( ...
                    psi, m_a, g_a, N_a, L_z, Mat_V, geom, hbar);

                vec_t(i_sampling) = i * dt;

                fig = Plot_single_component_state_sphere( ...
                    psi, geom.L_band, geom.method, ...
                    'Visible', 'off', ...
                    'FigureTitle', sprintf('$t = %.3e\\ \\mathrm{s}$', i*dt));

                fig_name = fullfile(temp_img_folder, sprintf('%05d.jpg', i_sampling));
                saveas(fig, fig_name);
                close(fig);

                i_sampling = i_sampling + 1;
            end
        end

        %--------------------------------------------------------------%
        % Trim unused preallocated entries
        %--------------------------------------------------------------%
        vec_mu  = vec_mu(1:i_sampling-1);
        vec_Ene = vec_Ene(1:i_sampling-1);
        vec_t   = vec_t(1:i_sampling-1);

        %--------------------------------------------------------------%
        % Final state figure
        %--------------------------------------------------------------%
        fig_final = Plot_single_component_state_sphere( ...
            psi, geom.L_band, geom.method, ...
            'Visible', 'off', ...
            'FigureTitle', sprintf('Prepared state at $t = %.3e\\ \\mathrm{s}$', n_iterations*dt));

        saveas(fig_final, fullfile(output_folder, 'final_state.fig'));
        saveas(fig_final, fullfile(output_folder, 'final_state.png'));
        close(fig_final);

        %--------------------------------------------------------------%
        % Diagnostics figure
        %--------------------------------------------------------------%
        fig_diag = figure('Visible', 'off', 'Position', [100, 100, 1200, 700]);

        tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        plot(vec_t, vec_mu, 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$\mu$ [J]', 'Interpreter', 'latex', 'FontSize', 15);
        title('$\mu(t)$', 'Interpreter', 'latex', 'FontSize', 16);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

        nexttile;
        plot(vec_t, vec_Ene, 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$E$ [J]', 'Interpreter', 'latex', 'FontSize', 15);
        title('$E(t)$', 'Interpreter', 'latex', 'FontSize', 16);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

        sgtitle('Imaginary-time diagnostics', 'Interpreter', 'latex', 'FontSize', 18);

        saveas(fig_diag, fullfile(output_folder, 'diagnostics.fig'));
        saveas(fig_diag, fullfile(output_folder, 'diagnostics.png'));
        close(fig_diag);

        %--------------------------------------------------------------%
        % Save output
        %--------------------------------------------------------------%
        output_mat_file = fullfile(output_folder, 'Single_component_monopole_on_sphere.mat');

        save(output_mat_file, ...
            'hbar', 'amu', 'a0', ...
            'm_a', 'N_a', ...
            'R_sphere', 'L_z', ...
            'a_aa', 'g_a', ...
            'g2D_est', 'n2D_est', 'mu_est', ...
            'n_iterations', 'dt', 'L_band', ...
            'geom', ...
            'q', ...
            'V0_barrier', 'sigma_barrier', ...
            'V0_pin', 'sigma_pin', ...
            'Mat_V', ...
            'psi', ...
            'vec_mu', 'vec_Ene', 'vec_t', ...
            'sample_frequency', 'params');

        fprintf('Saved MAT file: %s\n', output_mat_file);
        fprintf('Prepare_single_component_initial_state_sphere completed successfully in %.2f s.\n', toc(t_start));

        diary off;

    catch ME
        fprintf('\nERROR in Prepare_single_component_initial_state_sphere\n');
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        diary off;
        rethrow(ME);
    end
end


% ======================================================================= %
% Local helpers
% ======================================================================= %

function params = set_default_params(params)

    if ~isfield(params, 'N_a'),                  params.N_a = 5e5; end
    if ~isfield(params, 'R_sphere'),             params.R_sphere = 20e-6; end
    if ~isfield(params, 'L_z'),                  params.L_z = 2e-6 / sqrt(2*pi); end
    if ~isfield(params, 'a_aa'),                 params.a_aa = 52.0; end  % in units of a0

    if ~isfield(params, 'n_iterations_imag'),    params.n_iterations_imag = 1e3; end
    if ~isfield(params, 'dt_imag'),              params.dt_imag = 1e-6; end
    if ~isfield(params, 'L_band'),               params.L_band = 224; end
    if ~isfield(params, 'sample_frequency'),     params.sample_frequency = 1000; end

    if ~isfield(params, 'sigma_barrier'),        params.sigma_barrier = 0.05; end
    if ~isfield(params, 'V0_over_mu'),           params.V0_over_mu = 1.2; end

    % polar pinning potentials during imaginary-time preparation
    if ~isfield(params, 'use_pinning_poles'),    params.use_pinning_poles = true; end
    if ~isfield(params, 'V0_pin_over_mu'),       params.V0_pin_over_mu = 0.6; end
    if ~isfield(params, 'sigma_pin'),            params.sigma_pin = 0.05; end

    % initial-state smoothing / density shaping
    if ~isfield(params, 'delta_phase'),          params.delta_phase = 0.03; end
    if ~isfield(params, 'equator_notch_depth'),  params.equator_notch_depth = 0.95; end
    if ~isfield(params, 'seed_noise_amplitude'), params.seed_noise_amplitude = 0.0; end
end


function Mat_V = build_preparation_potential_sphere( ...
    geom, V0_barrier, sigma_barrier, V0_pin, sigma_pin)
% Preparation potential:
% - Gaussian barrier at the equator
% - Two Gaussian repulsive pinning potentials at North and South poles

    theta = geom.theta;

    % Equatorial barrier
    V_eq = V0_barrier * exp(-0.5 * ((theta - pi/2) ./ sigma_barrier).^2);

    % Polar pinning potentials
    V_N = V0_pin * exp(-0.5 * (theta ./ sigma_pin).^2);
    V_S = V0_pin * exp(-0.5 * ((theta - pi) ./ sigma_pin).^2);

    Mat_V = V_eq + V_N + V_S;
end


function psi = build_initial_patchwise_state_sphere(geom, q, params)
% Seed state with:
% - density approximately uniform in each hemisphere
% - strong density suppression near the equator
% - smooth interpolation of the phase sign near the equator

    theta = geom.theta;
    phi   = geom.phi;

    % Smooth sign switch across equator:
    % +1 in North, -1 in South
    sgn_smooth = tanh((pi/2 - theta) ./ params.delta_phase);

    % Patchwise phase:
    % North  -> +q*phi
    % South  -> -q*phi
    phase = q * sgn_smooth .* phi;

    % Density notch around the equator
    notch = params.equator_notch_depth * ...
        exp(-0.5 * ((theta - pi/2) ./ params.sigma_barrier).^2);

    rho0 = max(1e-12, 1 - notch);

    % Optional tiny random seed to break exact numerical symmetries
    if params.seed_noise_amplitude > 0
        noise = 1 + params.seed_noise_amplitude * ...
            (rand(size(theta)) - 0.5 + 1i*(rand(size(theta)) - 0.5));
    else
        noise = 1;
    end

    psi = sqrt(rho0) .* exp(1i * phase) .* noise;
end