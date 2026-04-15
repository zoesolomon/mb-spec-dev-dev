Infrared Spectroscopy with MB-pol

Codes included in this folder:

dipole: Calculates dipole moment of box using either 1B, 1B+NB, or 1B+2B+NB level (specify at runtime with either 1b, nb, or mb option)
   run mpirun -np 4 ./dipole <mb> <imcon> <dump frequency> <mu> ['{box info}'] <dump file>
	--- set mb to 1b for 1B monomer dipole moment surface only.
	--- set mb to nb for 1B DPS with included induction effects (1B+NB).
	--- set mb to mb for MB-mu DPS, consisting of explicit 1B and 2B effects and approximate NB effects (1B+2B+NB) (only for water).
        --- set mu to 0 to read the dipoles from a .lammpstrj, or 1 to calculate the 1B dipole from MB-mu DMS (only for water) 
	*** can be run in parallel over N cores with mpirun -np N option.
        *** dump file should be in a .lammpstrj format with each line as follows or in an .xyz format:
               ITEM: ATOMS id mol type q x y z v_dipx_perm v_dipy_perm v_dipz_perm v_dipnorm_perm v_dipx_ind v_dipy_ind v_dipz_ind v_dipnorm_ind
        *** if dump file is in an .xyz format, box info must be specified in the form '{xmin,xmax,ymin,ymax,zmin,zmax}'


response: Calculates average response function from one or multiple dipole moment data files
   ./response int_max dipole_surface(s) > response.dat
        --- set int_max to the maximum correlation time for your system.

spectrum: Calculates infrared spectrum from response function using cosine Fourier transform
   ./spectrum T volume tau alpha int_max < response.dat > spectrum.dat
        --- set T to be the temperature of the system.
        --- set tau to be the baseline shift.
        --- set alpha to between 0 and 1 for how much smoothing.
