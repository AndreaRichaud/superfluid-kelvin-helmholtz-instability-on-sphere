function oslo_custom = Load_oslo_colormap()
%LOAD_OSLO_COLORMAP Load the bundled Oslo density colormap.
%
% The repository includes the Oslo RGB table under
% src/visualization/colormaps/oslo.txt.
%
% OUTPUT
%   oslo_custom - 256-by-3 RGB colormap interpolating from white to the
%                 dark half of the inverted Oslo colormap.

utilities_directory = fileparts(mfilename('fullpath'));
src_directory = fileparts(utilities_directory);

oslo_file = fullfile( ...
    src_directory, ...
    'visualization', ...
    'colormaps', ...
    'oslo.txt');

if ~isfile(oslo_file)
    error('Bundled Oslo colormap file not found:\n  %s', oslo_file);
end

oslo = readmatrix(oslo_file);

if size(oslo, 2) ~= 3
    error(['Invalid Oslo colormap file: expected an N-by-3 RGB table, ' ...
        'but found a %d-by-%d array.'], ...
        size(oslo, 1), size(oslo, 2));
end

oslo = flipud(oslo);

n_colours = size(oslo, 1);
oslo_half = oslo(1:round(n_colours / 2), :);

white = [1, 1, 1];
target_colour = oslo_half(end, :);

s = linspace(0, 1, 256)';
oslo_custom = (1 - s) .* white + s .* target_colour;
end
