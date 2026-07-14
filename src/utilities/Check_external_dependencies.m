function Check_external_dependencies()
%CHECK_EXTERNAL_DEPENDENCIES Verify external MATLAB dependencies.
%
% The simulation requires the external SSHT MATLAB library.
% The Oslo RGB table used for visualization is bundled with this repository.

required_ssht_functions = {
    'ssht_sampling'
    'ssht_forward'
    'ssht_inverse'
    'ssht_plot_sphere'
};

for k = 1:numel(required_ssht_functions)
    if exist(required_ssht_functions{k}, 'file') == 0
        error(['Required SSHT function not found: %s.\n' ...
            'Install SSHT and add its MATLAB directory to the path.'], ...
            required_ssht_functions{k});
    end
end

% Sanity-check the bundled Oslo colormap before starting a long simulation.
Load_oslo_colormap;
end
