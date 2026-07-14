%% Run one vortex-charge case with oslo rendering
% Run the single-component real-time dynamics for one selected vortex
% charge. The corresponding prepared state must already exist in the case
% output directory.
%
% This script is a convenience entry point for regenerating a selected
% trajectory with custom real-time parameters. It is not a post-processing
% script.

clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));

Check_external_dependencies;

q = 8;

params = struct();
params.T_slope               = 0.05;
params.dt_real               = 1e-6;
params.tmax_real             = 0.10;
params.sample_frequency_real = 1000;
params.l_filter_real         = 18;

base_output_folder = fullfile(repository_root, 'output', 'Sweep_q');
output_folder = fullfile(base_output_folder, sprintf('q_%d', q));

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

Run_single_component_dynamics_sphere(q, params, output_folder);
