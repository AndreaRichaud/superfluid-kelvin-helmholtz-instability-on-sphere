function psi = sphere_normalize(psi, W)

    norm_psi = sqrt(sum(abs(psi(:)).^2 .* W(:)));
    psi = psi / norm_psi;

end