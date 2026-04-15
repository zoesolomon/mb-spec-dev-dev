Velocity-Velocity Power Spectroscopy

Codes included in this folder:

response: Calculates average response function from a trajectory file
   ./response <dump_file> <mb-spec.json> > response.dat
        --- set int_max to the maximum correlation time for your system.
        *** dump file should be in a .lammpstrj format with each line as follows or in an .xyz format:
               ITEM: ATOMS id mol type q vx vy vz 
        *** if dump file is in an .xyz format, box info must be specified in the form '{xmin,xmax,ymin,ymax,zmin,zmax}' 
        *** can be run in parallel over N cores with mpirun -np N option.       

spectrum: Calculates infrared spectrum from response function using cosine Fourier transform
   ./spectrum <response_file> <mb-spec.json> > spectrum.dat
        --- set T to be the temperature of the system.
        --- set tau to be the baseline shift.
        --- set alpha to between 0 and 1 for how much smoothing.
