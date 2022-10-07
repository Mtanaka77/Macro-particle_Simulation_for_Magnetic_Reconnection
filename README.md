Low-Frequency Electromagnetic Particle Simulation for High-Temperature Plasmas

Why is a current sheet of the distant magnetic tail of the earth become thinned and is a large energy released which is typically observed as magnetic reconnection ? There were many theories for the reconnection, in which Dr. Speicer had paid attention to the role of inertia resistivity of thinning the neutral sheet. Some time later, a particle simulation was executed to have shown the mechanism of earth's magnetotail reconnection which results in large energy release, in J. Geophysical Research (1995).

An electromagnetic particle simulation code is utilized for solar and magnetospheric space physics. Both electric and magnetic fields at low frequencies are solved by a slightly backward time decentering technique. Magnetic reconnection and the solar wind-earth magnetic field coupling are quite suitable for applying this simulation code.

One uses here the time decentered scheme in aimpl=0.6, while the time centered scheme in the explicit code (aimpl=0.5) is used in other directory of molecular dynamics simulations. It is noted, however, that finite errors in the divergence term accumulate which must be corrected if the finite difference coordinate space are utilized. Four physical units are, i) time: 1/wpe (c/wpe: electron inertia length), ii) length: c/wpe, iii) mass: electron mass, and iv) charge: electron charge. The program is written in Fortran 2003 and MPI version 3 for parallelization.

By the implicit scheme it is free from the Courant condition, that is, Dx(length)/Dt(time step) >< c, the speed of light. For the backward differential scheme in aimpl > 0.5, a time step may be dt~1.2/wpe to dump out plasma oscillations, while 2 \pi/dt wce > 1 for the kinetic ions and electrons. A large time step of 2 \pi/dt wce >> 1 is a good target of the drift-kinetic simulation, where a typical time step may be dt= 10/wpe.

One can enjoy simulations by changing system sizes and boundary conditions. For the present case, an equilibration of the pair of flux bundles is first tested in three dimensions. Kinetic ions and electrons are simulated in the igc=1 case specified in the rec_3d15A file. Then, one can start a merging of flux bundles. On the other case, the drift-kenetic electrons and kinetic ions are simulated as a large time step is used in the igc=2 case. But, one should note that heavy ions move kinetically while light electrons lose some of their particle freedom in the coordinate space.

Easy in-house plots are incorporated to check the current run in the simulation. Reading papers of this implicit particle simulation code (Ref. 1) and applications to magnetospheric space plasmas (Ref. 3) are highly recommended.


References:

1. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of three dimensions, J.Comput. Phys., 107, 124-145 (1993).

2. M. Tanaka, Macro-EM particle simulation method and a study of collisionless magnetic reconnection, Comput.Phys.Commun., 87, 117-138 (1995).

3. M. Tanaka, Macro-particle simulations of collisionless magnetic reconnection, Phys.Plasmas, 2, 2920-2930 (1995).

4. H. Shimazu, M. Tanaka, and S. Machida, The behavior of heavy ions in collisionless parallel shocks generated by the solar wind and planetary plasma interactions, J.Geophys.Res., 101, 27565-27571 (1996).


