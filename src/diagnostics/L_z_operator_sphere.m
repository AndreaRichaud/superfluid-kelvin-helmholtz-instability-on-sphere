function Lz_psi = L_z_operator_sphere(psi, geom, hbar)
%L_Z_OPERATOR_SPHERE
% Apply the angular momentum operator L_z = -i*hbar*d/dphi
% to a wavefunction on the sphere.

    reality_flag = isreal(psi);

    flm = ssht_forward(psi, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality_flag);

    flm_Lz = hbar * geom.m_vec(:) .* flm(:);

    Lz_psi = ssht_inverse(flm_Lz, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', false);
end