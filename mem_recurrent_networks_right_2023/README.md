# Memory of recurrent networks: Do we compute it right? (2023)

### Giovanni Ballarin<sup>1</sup>, Lyudmila Grigoryeva<sup>2,3</sup>, Juan-Pablo Ortega<sup>4</sup>

1. Department of Economics, University of Mannheim
2. Faculty of Mathematics and Statistics, University of St. Gallen
3. Honorary Associate Professor, Department of Statistics, University of Warwick
4. Division of Mathematical Sciences, Nanyang Technological University

Links to paper: [[ArXiv]](https://arxiv.org/abs/2305.01457) [[ResearchGate]](https://www.researchgate.net/publication/370462485_Memory_of_recurrent_networks_Do_we_compute_it_right)

Main scripts:

- `mc_sample_floor_effect.m`: reproduces Figure 1;

- `mc_gamma_eigs_plot.m`: reproduces plots in Figure 2;

- `mc_krylov_eigs.m`: reproduces plots in Figure 3 by choosing the appropriate distributions for `A`;

- `mc_compute_bands.m`: reproduces plots in Figures 4 and 5 by choosing the appropriate distributions for `A` and `C`;

- `mc_plot_connect_eigenvalues.m`: reproduces plots in Figure 6 by choosing the appropriate distribution for `A`;

- `mc_compute_showcase.m`: allows to create simple plots for comparing the naive MC estimation method with OSM and OSM+, our proposals;
