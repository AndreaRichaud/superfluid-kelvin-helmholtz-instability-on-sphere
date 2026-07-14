%% Build comparison videos from Oslo frames for q_1 ... q_9
clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));

base_folder = fullfile(repository_root, 'output', 'Sweep_q');

q_values = 1:9;

target_duration_sec = 10;
video_quality = 100;

%% Reference time vector
ref_file = fullfile(base_folder, 'q_1', ...
    'Real_time_dynamics_single_component_on_sphere_oslo.mat');

if ~exist(ref_file, 'file')
    error('Reference MAT file not found:\n%s', ref_file);
end

S = load(ref_file, 'vec_t');
vec_t = S.vec_t(:);

%% Load all frame lists from Temp_real_oslo
all_frames = cell(numel(q_values),1);
min_frames = inf;

for k = 1:numel(q_values)

    q_folder = fullfile(base_folder, sprintf('q_%d', q_values(k)));
    temp_folder = fullfile(q_folder, 'Temp_real_oslo');

    if ~exist(temp_folder, 'dir')
        error('Missing Temp_real_oslo in %s', q_folder);
    end

    frame_files = dir(fullfile(temp_folder, '*.jpg'));

    if isempty(frame_files)
        error('No jpg frames in %s', temp_folder);
    end

    [~, idx] = sort({frame_files.name});
    frame_files = frame_files(idx);

    all_frames{k} = frame_files;
    min_frames = min(min_frames, numel(frame_files));
end

min_frames = min(min_frames, numel(vec_t));

fprintf('Using %d frames.\n', min_frames);

%% =========================================================
% Density video: left half of each frame
%% =========================================================

make_comparison_video( ...
    base_folder, ...
    q_values, ...
    all_frames, ...
    vec_t, ...
    min_frames, ...
    target_duration_sec, ...
    video_quality, ...
    'left', ...
    fullfile(base_folder, 'comparison_q1_to_q9_density_oslo.mp4'));

%% =========================================================
% Phase video: right half of each frame
%% =========================================================

make_comparison_video( ...
    base_folder, ...
    q_values, ...
    all_frames, ...
    vec_t, ...
    min_frames, ...
    target_duration_sec, ...
    video_quality, ...
    'right', ...
    fullfile(base_folder, 'comparison_q1_to_q9_phase_oslo.mp4'));

%% =========================================================
% Local function
%% =========================================================

function make_comparison_video(base_folder, q_values, all_frames, vec_t, ...
    min_frames, target_duration_sec, video_quality, panel_side, output_video_file)

    fps = max(1, round(min_frames / target_duration_sec));

    v = VideoWriter(output_video_file, 'MPEG-4');
    v.FrameRate = fps;
    v.Quality = video_quality;
    open(v);

    for j = 1:min_frames

        imgs = cell(numel(q_values),1);

        for k = 1:numel(q_values)

            q_folder = fullfile(base_folder, sprintf('q_%d', q_values(k)));
            temp_folder = fullfile(q_folder, 'Temp_real_oslo');

            img_path = fullfile(temp_folder, all_frames{k}(j).name);
            img = imread(img_path);

            width = size(img,2);

            switch lower(panel_side)

                case 'left'
                    img = img(:, 1:floor(width/2), :);

                case 'right'
                    img = img(:, floor(width/2)+1:end, :);

                otherwise
                    error('panel_side must be either left or right.')

            end

            % Remove local t label / sgtitle region
            h_mask = round(0.10 * size(img,1));
            img(1:h_mask, :, :) = 255;

            imgs{k} = img;
        end

        % Resize all panels to the same size
        target_size = size(imgs{1});

        for k = 1:numel(q_values)
            imgs{k} = imresize(imgs{k}, [target_size(1), target_size(2)]);
        end

        % 3x3 grid
        row1 = [imgs{1}, imgs{2}, imgs{3}];
        row2 = [imgs{4}, imgs{5}, imgs{6}];
        row3 = [imgs{7}, imgs{8}, imgs{9}];

        grid_img = [row1; row2; row3];

        % Global time label
        t_current = vec_t(j);

        grid_img = insertText(grid_img, [20 20], ...
            sprintf('t = %.4f s', t_current), ...
            'FontSize', 40, ...
            'BoxColor', 'white', ...
            'TextColor', 'black');

        writeVideo(v, grid_img);

        if mod(j,100) == 0
            fprintf('%s: frame %d / %d\n', output_video_file, j, min_frames);
        end
    end

    close(v);

    fprintf('Saved video:\n%s\n', output_video_file);

end