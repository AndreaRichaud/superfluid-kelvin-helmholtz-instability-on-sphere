function psi_filt = Spectral_filter_sphere(psi, geom, p)

    reality_flag = isreal(psi);

    flm = ssht_forward(psi, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality_flag);

    sigma = zeros(size(flm));

    ind = 1;
    for ell = 0:(geom.L_band-1)
        filt = exp(-(ell/(geom.L_band-1))^p);
        for m = -ell:ell
            sigma(ind) = filt;
            ind = ind + 1;
        end
    end

    flm = sigma .* flm;

    psi_filt = ssht_inverse(flm, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality_flag);
end