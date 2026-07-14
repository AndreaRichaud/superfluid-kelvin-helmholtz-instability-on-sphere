function lap_f = spherical_laplacian(f, geom, varargin)
%SPHERICAL_LAPLACIAN Laplace-Beltrami operator on a sphere
%
%   lap_f = spherical_laplacian(f, geom)
%   lap_f = spherical_laplacian(f, geom, 'Reality', true)
%
% INPUTS
%   f       : function sampled on the SSHT grid stored in geom
%   geom    : geometry struct produced by Build_geometry_sphere, containing:
%               geom.R_sphere
%               geom.L_band
%               geom.method
%               geom.lap_eigs
%
% OPTIONAL NAME-VALUE PAIRS
%   'Reality' : true if f is real-valued (default: isreal(f))
%
% OUTPUT
%   lap_f   : Laplace-Beltrami of f on a sphere of radius geom.R_sphere
%
% NOTES
%   The operator implemented is
%
%       Delta_surf = (1/R^2) * Delta_{S^2}
%
%   and in spherical-harmonic space:
%
%       Delta_{S^2} Y_{lm} = -l(l+1) Y_{lm}.
%
%   Therefore:
%
%       Delta_surf Y_{lm} = -l(l+1)/R^2 * Y_{lm}.
%

    p = inputParser;
    addRequired(p, 'f', @(x) isnumeric(x));
    addRequired(p, 'geom', @(x) isstruct(x) && ...
        isfield(x,'L_band') && ...
        isfield(x,'method') && ...
        isfield(x,'lap_eigs'));
    addParameter(p, 'Reality', isreal(f), @(x) islogical(x) || isnumeric(x));
    parse(p, f, geom, varargin{:});

    reality = logical(p.Results.Reality);

    % Basic consistency check
    if ~isequal(size(f), size(geom.theta))
        error('spherical_laplacian:SizeMismatch', ...
            'Input f must have the same size as geom.theta / geom.phi.');
    end

    % Forward spherical harmonic transform
    flm = ssht_forward(f, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality);

    % Apply spectral Laplacian: -l(l+1)/R^2
    % geom.lap_eigs stores +l(l+1)/R^2, so we multiply by -geom.lap_eigs
    flm_lap = -geom.lap_eigs(:) .* flm(:);

    % Inverse transform
    lap_f = ssht_inverse(flm_lap, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality);
end