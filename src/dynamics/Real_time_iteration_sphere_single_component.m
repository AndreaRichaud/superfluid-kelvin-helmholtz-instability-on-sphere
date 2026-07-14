function Hpsi = Real_time_iteration_sphere_single_component( ...
    psi, ...
    m_a, ...
    g_a, ...
    N_a, ...
    L_z, ...
    Mat_V, ...
    geom, hbar)

% REAL_TIME_ITERATION_SPHERE_SINGLE_COMPONENT
% Compute H psi for a single-component real-time GPE on a sphere.
%
% INPUTS
%   psi      : wavefunction
%   m_a      : mass
%   g_a      : intra-species 3D coupling
%   N_a      : atom number
%   L_z      : effective radial thickness used for 3D -> 2D reduction
%   Mat_V    : external potential
%   geom     : spherical geometry struct
%   hbar     : reduced Planck constant
%
% OUTPUT
%   Hpsi     : Hamiltonian action on psi

    % Effective 2D coupling
    g_a_2D = g_a / L_z;

    % Laplace-Beltrami term
    Laplacian = spherical_laplacian(psi, geom, 'Reality', isreal(psi));

    % Nonlinear intra-species term
    Intra = g_a_2D * N_a * abs(psi).^2 .* psi;

    % External potential
    Pot = Mat_V .* psi;

    % Hamiltonian action
    Hpsi = -(hbar^2/(2*m_a)) * Laplacian + Intra + Pot;

end