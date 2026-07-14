function Run_single_component_dynamics_sphere(q, params, output_folder)
%RUN_SINGLE_COMPONENT_DYNAMICS_SPHERE Evolve the prepared state in real time.
%
% Load the prepared single-component state, lower the equatorial barrier
% linearly over T_slope, and evolve the Gross-Pitaevskii equation on the
% sphere. The routine saves diagnostic data, JPEG frames rendered with the
% oslo colormap, and wavefunction snapshots used by the vortex-tracking
% post-processing workflow.
%
% INPUTS
%   q             - integer vortex charge, used for validation and labels
%   params        - real-time parameter structure
%   output_folder - case directory containing the prepared-state MAT file


    if nargin < 2 || isempty(params)
        params = struct();
    end

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    temp_real_folder = fullfile(output_folder, 'Temp_real_oslo');
    if ~exist(temp_real_folder, 'dir')
        mkdir(temp_real_folder);
    end

    psi_save_folder = fullfile(output_folder, 'psi_a_vs_t');
    if ~exist(psi_save_folder, 'dir')
        mkdir(psi_save_folder);
    end

    log_file = fullfile(output_folder, 'log_real_time_dynamics_single_component.txt');
    if exist(log_file, 'file')
        delete(log_file);
    end

    diary(log_file);
    diary on;

    t_start = tic;

    try
        fprintf('============================================================\n');
        fprintf('Run_single_component_dynamics_sphere started\n');
        fprintf('q = %d\n', q);
        fprintf('Output folder: %s\n', output_folder);
        fprintf('Start time: %s\n', datestr(now));
        fprintf('============================================================\n');

        %% Path setup
        this_file_folder = fileparts(mfilename('fullpath'));
        addpath(genpath(this_file_folder));

        %% Load prepared initial state
        input_mat_file = fullfile(output_folder, ...
            'Single_component_monopole_on_sphere.mat');

        if ~exist(input_mat_file, 'file')
            error('Input file not found: %s', input_mat_file);
        end

        S = load(input_mat_file);

        if isfield(S, 'params')
            params = inherit_missing_fields(params, S.params);
        end

        params = set_default_movie_params(params);

        hbar          = S.hbar;
        m_a           = S.m_a;
        N_a           = S.N_a;
        L_z           = S.L_z;
        g_a           = S.g_a;
        geom          = S.geom;
        psi           = S.psi;
        V0_barrier    = S.V0_barrier;
        sigma_barrier = S.sigma_barrier;

        if isfield(S, 'q') && S.q ~= q
            warning(['Loaded q (%d) differs from function input q (%d). ' ...
                     'Using function input q only for labeling.'], S.q, q);
        end

        fprintf('Prepared-state file loaded successfully.\n');

        %% Real-time evolution parameters
        dt               = params.dt_real;
        tmax             = params.tmax_real;
        T_slope          = params.T_slope;
        sample_frequency = params.sample_frequency_real;
        l_filter         = params.l_filter_real;

        n_iterations = round(tmax / dt);
        n_samples    = ceil(n_iterations / sample_frequency) + 1;

        vec_t    = zeros(n_samples,1);
        vec_Ene  = zeros(n_samples,1);
        vec_norm = zeros(n_samples,1);
        vec_V0   = zeros(n_samples,1);
        vec_L_z  = zeros(n_samples,1);

        c_p = 1;

        fprintf('Real-time evolution setup completed.\n');
        fprintf('n_iterations = %d\n', n_iterations);
        fprintf('dt_real      = %.3e s\n', dt);
        fprintf('tmax_real    = %.3e s\n', tmax);
        fprintf('T_slope      = %.3e s\n', T_slope);
        fprintf('sample_frequency_real = %d\n', sample_frequency);

        %% Real-time evolution loop
        for i = 1:n_iterations

            t_now = (i-1) * dt;

            if t_now <= T_slope
                V0_t = V0_barrier * (1 - t_now / T_slope);
            else
                V0_t = 0;
            end

            Mat_V = build_equatorial_barrier_sphere( ...
                geom, V0_t, sigma_barrier);

            psi = Real_time_evolve_sphere_single_component( ...
                psi, ...
                m_a, ...
                g_a, ...
                N_a, ...
                L_z, ...
                Mat_V, ...
                geom, hbar, dt, l_filter);

            if i == 1 || mod(i, sample_frequency) == 0

                elapsed_time = toc(t_start);
                t_sample = i * dt;

                vec_t(c_p)    = t_sample;
                vec_V0(c_p)   = V0_t;
                vec_norm(c_p) = real(sum(abs(psi(:)).^2 .* geom.W(:)));

                vec_Ene(c_p) = Compute_total_energy_single_component_sphere( ...
                    psi, ...
                    m_a, ...
                    g_a, ...
                    N_a, ...
                    L_z, ...
                    Mat_V, ...
                    geom, hbar);

                vec_L_z(c_p) = Compute_total_angular_momentum_z_single_component_sphere( ...
                    psi, N_a, geom, hbar);

                if c_p == 1

                    dE_percent  = 0;
                    dLz_percent = 0;

                else

                    E0  = vec_Ene(1);
                    Lz0 = vec_L_z(1);

                    if abs(E0) > eps
                        dE_percent = 100 * (vec_Ene(c_p) - E0) / E0;
                    else
                        dE_percent = NaN;
                    end

                    if abs(Lz0) > eps
                        dLz_percent = 100 * (vec_L_z(c_p) - Lz0) / Lz0;
                    else
                        dLz_percent = NaN;
                    end

                end

                fprintf(['Progress = %.2f %% in %.0f s | ' ...
                    'V0(t) = %.3e J | norm = %.8f | ' ...
                    'dE = %+10.3e %% | dL_z = %+10.3e %%\n'], ...
                    100*i/n_iterations, elapsed_time, vec_V0(c_p), vec_norm(c_p), ...
                    dE_percent, dLz_percent);

                %% Oslo frame
                fig = Plot_single_component_state_oslo_sphere( ...
                    psi, geom.L_band, geom.method, ...
                    'Visible', 'off', ...
                    'FigureTitle', sprintf('$t = %.4f\\ \\mathrm{s}$', t_sample));

                fig_name = fullfile(temp_real_folder, sprintf('%05d.jpg', c_p));
                saveas(fig, fig_name);
                close(fig);

                %% Save wavefunction snapshot
                psi_a = psi; %#ok<NASGU>
                t = t_sample; %#ok<NASGU>

                psi_file = fullfile(psi_save_folder, ...
                    sprintf('psi_a_t_%05d.mat', c_p));

                save(psi_file, ...
                    'psi_a', ...
                    't', ...
                    'q', ...
                    'geom', ...
                    'hbar', ...
                    'm_a', ...
                    'N_a', ...
                    'L_z', ...
                    'g_a');

                c_p = c_p + 1;

            end
        end

        %% Trim arrays
        c_p = c_p - 1;

        vec_t    = vec_t(1:c_p);
        vec_Ene  = vec_Ene(1:c_p);
        vec_norm = vec_norm(1:c_p);
        vec_V0   = vec_V0(1:c_p);
        vec_L_z  = vec_L_z(1:c_p);

        vec_t    = vec_t(:);
        vec_Ene  = vec_Ene(:);
        vec_norm = vec_norm(:);
        vec_V0   = vec_V0(:);
        vec_L_z  = vec_L_z(:);

        V0_final = max(V0_barrier * (1 - n_iterations*dt / T_slope), 0);

        Mat_V_final = build_equatorial_barrier_sphere( ...
            geom, V0_final, sigma_barrier);

        %% Save real-time dynamics MAT
        output_mat_file = fullfile(output_folder, ...
            'Real_time_dynamics_single_component_on_sphere_oslo.mat');

        save(output_mat_file, ...
            'q', ...
            'n_iterations', 'sample_frequency', ...
            'dt', 'tmax', 'T_slope', ...
            'psi', ...
            'V0_barrier', 'V0_final', 'sigma_barrier', ...
            'Mat_V_final', ...
            'vec_t', 'vec_Ene', 'vec_norm', 'vec_V0', 'vec_L_z', ...
            'm_a', 'N_a', 'L_z', ...
            'g_a', ...
            'geom', 'hbar', ...
            'params');

        fprintf('Saved MAT file: %s\n', output_mat_file);

        %% Final state figure
        fig_final = Plot_single_component_state_oslo_sphere( ...
            psi, geom.L_band, geom.method, ...
            'Visible', 'off', ...
            'FigureTitle', sprintf('Final state at $t = %.4f\\ \\mathrm{s}$', ...
            n_iterations*dt));

        saveas(fig_final, fullfile(output_folder, 'final_real_time_state_oslo.fig'));
        saveas(fig_final, fullfile(output_folder, 'final_real_time_state_oslo.png'));
        close(fig_final);

        %% Diagnostics figure
        fig_diag = figure('Visible','off','Position',[100 100 900 1100]);

        tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

        nexttile;
        plot(vec_t, vec_Ene, 'k', 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
        ylabel('$E(t)$ [J]','Interpreter','latex','FontSize',16);
        title('$E(t)$','Interpreter','latex','FontSize',18);
        set(gca,'TickLabelInterpreter','latex','FontSize',14);

        nexttile;
        plot(vec_t, vec_L_z, 'k', 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
        ylabel('$L_z(t)$ [J\,s]','Interpreter','latex','FontSize',16);
        title('$L_z(t)$','Interpreter','latex','FontSize',18);
        set(gca,'TickLabelInterpreter','latex','FontSize',14);

        nexttile;
        plot(vec_t, vec_norm, 'k', 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
        ylabel('$\int |\psi|^2 d\Omega$','Interpreter','latex','FontSize',16);
        title('Norm conservation','Interpreter','latex','FontSize',18);
        set(gca,'TickLabelInterpreter','latex','FontSize',14);

        nexttile;
        plot(vec_t, vec_V0, 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
        ylabel('$V_0(t)$ [J]','Interpreter','latex','FontSize',16);
        title('Barrier ramp','Interpreter','latex','FontSize',18);
        set(gca,'TickLabelInterpreter','latex','FontSize',14);

        sgtitle('Real-time diagnostics','Interpreter','latex','FontSize',18);

        saveas(fig_diag, fullfile(output_folder, 'real_time_diagnostics_oslo.fig'));
        saveas(fig_diag, fullfile(output_folder, 'real_time_diagnostics_oslo.png'));
        close(fig_diag);

        fprintf('Run_single_component_dynamics_sphere completed successfully in %.2f s.\n', toc(t_start));

        diary off;

    catch ME

        fprintf('\nERROR in Run_single_component_dynamics_sphere\n');
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        diary off;
        rethrow(ME);

    end
end

% ======================================================================= %
% Local helpers
% ======================================================================= %

function params = set_default_movie_params(params)

    if ~isfield(params, 'dt_real'),               params.dt_real = 1e-6; end
    if ~isfield(params, 'tmax_real'),             params.tmax_real = 0.10; end
    if ~isfield(params, 'T_slope'),               params.T_slope = 0.05; end
    if ~isfield(params, 'sample_frequency_real'), params.sample_frequency_real = 1000; end
    if ~isfield(params, 'l_filter_real'),         params.l_filter_real = 18; end

end

function params = inherit_missing_fields(params, params_from_file)

    f = fieldnames(params_from_file);

    for k = 1:numel(f)

        if ~isfield(params, f{k})
            params.(f{k}) = params_from_file.(f{k});
        end

    end

end

function Mat_V = build_equatorial_barrier_sphere(geom, V0, sigma_barrier)

    theta = geom.theta;
    Mat_V = V0 * exp(-0.5 * ((theta - pi/2) ./ sigma_barrier).^2);

end