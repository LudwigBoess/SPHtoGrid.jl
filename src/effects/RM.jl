function rotation_measure(n_cm3::Real, B_par::Real)
    return 812 * n_cm3 * 1.e3 * B_par * 1.e6
end