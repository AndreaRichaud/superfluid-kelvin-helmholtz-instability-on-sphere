function Ene = Compute_total_energy_single_component_sphere( ...
    psi, ...
    m_a, ...
    g_a, ...
    N_a, ...
    L_z, ...
    Mat_V, ...
    geom, hbar)

% ENERGY_COMPUTATOR_SPHERE_SINGLE_COMPONENT
% Compute the total energy of a single-component GPE on a sphere.
%
% INPUTS
%   psi      : wavefunction normalized to 1 on the sphere
%   m_a      : mass
%   g_a      : intra-species 3D coupling
%   N_a      : atom number
%   L_z      : effective radial thickness used for 3D -> 2D reduction
%   Mat_V    : external potential acting on the component
%   geom     : spherical geometry struct
%   hbar     : reduced Planck constant
%
% OUTPUT
%   Ene      : total energy [J]

    % Effective 2D coupling
    g_a_2D = g_a / L_z;

    % Laplace-Beltrami term
    Laplacian = spherical_laplacian(psi, geom, 'Reality', isreal(psi));

    % Kinetic energy
    Kinetic = -N_a * hbar^2/(2*m_a) * sum(conj(psi(:)) .* Laplacian(:) .* geom.W(:));

    % Potential energy
    Pot = N_a * sum(Mat_V(:) .* abs(psi(:)).^2 .* geom.W(:));

    % Intra-species interaction energy
    Intra = g_a_2D * N_a^2 / 2 * sum(abs(psi(:)).^4 .* geom.W(:));

    % Total energy
    Ene = real(Kinetic + Pot + Intra);

end