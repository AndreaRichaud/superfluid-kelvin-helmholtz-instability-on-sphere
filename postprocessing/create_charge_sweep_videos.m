%% Build videos from Temp_real_oslo folders inside Sweep_q
clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));

base_folder = fullfile(repository_root, 'output', 'Sweep_q');

% Find all q_* folders
dir_list = dir(fullfile(base_folder, 'q_*'));
dir_list = dir_list([dir_list.isdir]);

target_duration_sec = 10;
video_quality = 100;

for k = 1:numel(dir_list)

    case_folder = fullfile(base_folder, dir_list(k).name);
    temp_real_folder = fullfile(case_folder, 'Temp_real_oslo');

    if ~exist(temp_real_folder, 'dir')
        fprintf('Skipping %s: Temp_real_oslo not found\n', dir_list(k).name);
        continue
    end

    frame_files = dir(fullfile(temp_real_folder, '*.jpg'));

    if isempty(frame_files)
        fprintf('Skipping %s: no jpg frames found in Temp_real_oslo\n', dir_list(k).name);
        continue
    end

    [~, idx] = sort({frame_files.name});
    frame_files = frame_files(idx);

    n_frames = numel(frame_files);

    fps = max(1, round(n_frames / target_duration_sec));

    output_video_file = fullfile(case_folder, ...
        sprintf('%s_oslo.mp4', dir_list(k).name));

    fprintf('Creating Oslo video for %s\n', dir_list(k).name);
    fprintf('  Frames: %d\n', n_frames);
    fprintf('  FPS: %d\n', fps);
    fprintf('  Output: %s\n', output_video_file);

    v = VideoWriter(output_video_file, 'MPEG-4');
    v.FrameRate = fps;
    v.Quality = video_quality;
    open(v);

    for j = 1:n_frames

        img_path = fullfile(temp_real_folder, frame_files(j).name);
        img = imread(img_path);

        writeVideo(v, img);

    end

    close(v);

    actual_duration = n_frames / fps;
    fprintf('  Done. Video duration: %.2f s\n\n', actual_duration);

end