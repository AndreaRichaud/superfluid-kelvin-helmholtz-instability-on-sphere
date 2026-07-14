function Lz = Compute_total_angular_momentum_z_single_component_sphere( ...
    psi, N_a, geom, hbar)
%TOTAL_ANGULAR_MOMENTUM_L_Z_SPHERE_SINGLE_COMPONENT
% Compute the angular momentum along z for a single-component condensate
% on the sphere.
%
% INPUTS
%   psi      : wavefunction normalized to 1
%   N_a      : atom number
%   geom     : geometry struct
%   hbar     : reduced Planck constant
%
% OUTPUT
%   Lz       : total angular momentum along z

    % Apply angular momentum operator
    Lz_psi = L_z_operator_sphere(psi, geom, hbar);

    % Expectation value
    Lz = N_a * real(sum(conj(psi(:)) .* Lz_psi(:) .* geom.W(:)));

end