chapter_2 = ["stellar_minimal"]
chapter_3 = ["plot_plummer", "plot_hydro_density_distributions", "salpeter", "sun_venus_earth", "make_binary_star", ]
chapter_4 = ["gravity_minimal", "sun_venus_earth", "plot_gravity", "anim_gravity", "montgomery", "hierarchical_quadruple", ]
chapter_5 = ["calculate_core_temperature_density", "initialize_single_star", "stellar_minimal", "multiple_stellar_code", "merge_two_stars", "binary_evolution", "massloss_response", ]
chapter_6 = ["hydro_minimal", "hydro_outflow_particles", "hydro_sink_particles", "hydro_disk_with_bump", "explode_evolved_star", "helium_star", ]
chapter_7 = ["rad_minimal", "radiative_h2cloud", ]
chapter_8 = ["gravity_stellar_minimal", "gravity_stellar_collision", "plot_fractal_clumpy_cluster", "gravity_stellar_eventdriven", "merge_two_stars_and_evolve", "merge_two_stars_sph", "gravity_kepler_disks", "basic_multiples", "kira", "lonely_planets", ]
chapter_9 = ["gravity_potential", "gravity_drifter", "two_body_bridge", "three_body_bridge", "three_body_bridge_hierarchical", "gravity_hydro", "relax_gas_and_stars", "evolve_xi_tau_primary", ]
chapter_10 = ["multiples_minimal", "kira", ]  
appendix_a = ["multiple_stellar_threaded", ]
appendix_b = []
appendix_c = ["class_example", ]

examples_per_chapter = {
    "2": chapter_2,
    "fundamental": chapter_2,
    "Fundamental Concepts": chapter_2,
    "3": chapter_3,
    "ic": chapter_3,
    "Initial Conditions": chapter_3,
    "4": chapter_4,
    "Gravity": chapter_4,
    "Gravitational Dynamics": chapter_4,
    "5": chapter_5,
    "stellar": chapter_5,
    "Stellar Structure and Evolution": chapter_5,
    "6": chapter_6,
    "hydro": chapter_6,
    "Hydrodynamics": chapter_6,
    "7": chapter_7,
    "radiative": chapter_7,
    "Radiative Transfer": chapter_7,
    "8": chapter_8,
    "coupling": chapter_8,
    "Elementary Coupling Strategies": chapter_8,
    "9": chapter_9,
    "bridge": chapter_9,
    "Advanced Coupling Strategies": chapter_9,
    "10": chapter_10,
    "Expanding and extending AMUSE": chapter_10,
    "a": appendix_a,
    "AMUSE Fundamentals": appendix_a,
    "b": appendix_b,
    "How to add a code to AMUSE": appendix_b,
    "c": appendix_c,
    "Programming Primer": appendix_c,
}
