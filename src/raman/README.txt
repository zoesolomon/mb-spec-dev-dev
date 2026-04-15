Raman Spectroscopy with MB-pol

Codes included in this folder:

polarizability: Calculates polarizability tensor of box using either 1B, 1B+NB, or 1B+2B+NB level (specify at runtime with either 1b, nb, or mb option)
   mpirun -np N ./polarizability mb imcon dt dump.dipoles > pol.dat
	--- set mb to 1b for 1B Avila monomer polarizability tensor only.
	--- set mb to nb for 1B tensor with included induction effects (1B+NB). ***very expensive computationally
	--- set mb to mb for MB-alpha, consisting of explicit 1B and 2B effects and approximate NB effects (1B+2B+NB). ***very expensive computationally
	*** can be run in parallel over N cores with mpirun -np N option

response_isotropic: Calculates average isotropic response function from one or multiple polarizability tensor data files
   ./response_isotropic pol1.dat pol2.dat ... > resp.dat

response_depolarized: Calculates average depolarized response function from one or multiple polarizability tensor data files
   ./response_depolarized pol1.dat pol2.dat ... > resp.dat

spectrum: Calculates Raman spectrum from response function using cosine Fourier transform
   ./spectrum volume tau alpha max_time < resp.dat > spectrum.dat
