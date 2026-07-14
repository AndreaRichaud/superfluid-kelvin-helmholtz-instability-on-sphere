%% Sweep over the vortex charge q
% Run the single-component spherical Gross-Pitaevskii workflow for q = 1:10.
%
% By default, the script prepares each initial state from scratch before
% running the real-time dynamics. Set prepare_initial_states = false to
% reuse existing prepared-state files.
%
% Output is written below output/Sweep_q/q_<value>/.
%
% External requirements:
%   - SSHT MATLAB library
%   - Parallel Computing Toolbox

clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));

Check_external_dependencies;

q_values = 1:10;
prepare_initial_states = true;

params = struct();
params.T_slope               = 0.05;
params.dt_real               = 1e-6;
params.tmax_real             = 0.15;
params.sample_frequency_real = 1000;
params.l_filter_real         = 18;

output_root = fullfile(repository_root, 'output');
base_output_folder = fullfile(output_root, 'Sweep_q');

if ~exist(base_output_folder, 'dir')
    mkdir(base_output_folder);
end

parfor k = 1:numel(q_values)

    q = q_values(k);
    output_folder = fullfile(base_output_folder, sprintf('q_%d', q));

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    if prepare_initial_states
        Prepare_single_component_initial_state_sphere( ...
            q, params, output_folder);
    end

    Run_single_component_dynamics_sphere( ...
        q, params, output_folder);
end
