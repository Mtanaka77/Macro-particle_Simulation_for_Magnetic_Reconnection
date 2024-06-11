## Macro-Particle Simulation for Magnetic Reconnection ## 

This page is discussed on the largescale electromagnetic particle simulation code  
(J. Comp. Physics, Tanaka, 1993), magnetic reconnection (Phys. Plasmas, Tanaka, 1995)
and the solar wind and planetary plasma interactions in collisionless parallel shocks (J.Geophys.Res., 1996).
. 
### Magnetic Reconnection in Solar-Magnetospheric Couplings ###

Why is a large amount of the solar-eartth energy released intermittently 
in the distant magnetotail ? This energy release is suddenly and typically observed 
as magnetic reconnection. 

There were many theories for magnetic reconnection including from classical Dungey's 
theory to nuclear-fusion oriented anomalous resistivity. To our view, it was noted 
that Dr. Speicer said as 'hypothesis' inertia resistivity of thinning the current sheet,
which was not paid attention because of popularity of anomalous resistivity. 

Much later by using the particle-in-cell simulation, it was clearly proved 
that 'inertia of ions and electrons' is the key of input and output flows of 
magnetic reconnection resulting in the large energy release of earth's magnetotail (Ref.1).

### Implicit Particle-in-Cell Simulation Code ###

An electromagnetic particle simulation code is utilized in solar and magnetospheric 
space physics (Ref. 1, 4-5). 
Both electric and magnetic fields are solved at low frequencies by a slightly 
backward time decentering technique (Ref. 2,3). 
The backward de-centering code does not affect low-frequency phenomena, 
\omega*Dt << 1 with Dt the time step and \omega the inverse of electron plasma 
frequency (JCP, 1993).
Magnetic reconnection and the solar wind-earth's magnetic field coupling 
are quite suitable for applying this simulation code.

In the simulation, one utilizes the time-decentered scheme that is typically \aimpl=0.6, 
while the time-centered scheme in the explicit code (\aimpl=0.5) is used 
only in other directories of molecular dynamics simulations of \omega Dt <<1.

Four physical units in the sumulation are, i) time: 1/wpe, ii) length: c/wpe 
 (where c/wpe is the electron inertia length), iii) mass: electron mass, and 
iv) charge: electron charge. 
The program is written in Fortran 2003/2008 and it is coded for parallelization 
by MPI version 3 or 4.
The title, major references and remarks of this simulation code are written 
in the top part of the @mrg37_023A.f03 file.
The major subroutines are named as /fulmov/, /emfild/ and /cfpsol/, which are used 
in every time steps, while /escorr/ and /fulmv2/ are called in 5 time steps interval. 

The correction to the longitudinal part of the electric field is made in /escorr/. 
Although the Poisson's equation for the electric field is to be solved only initially, 
it is actually not quite true in the Maxwell equation since numerical errors 
do accumulate in time steps (see Ref. 2,3).

Supporting subroutines used are /partpc/ and /srimp1/ - /srimp4/. 
Important these subroutines are precisely explained as comments.
Two additional files are necessary with the code number  '23', 
as the parameter file param_A23A.h and the configure file rec_3d23A.

By the implicit scheme it is free from the Courant condition, that is, Dx(length)/Dt(time step) >< c, the speed of light. For the backward differential scheme in aimpl > 0.5, a time step may be dt~1.2/wpe in order to dump out plasma oscillations at plasma frequency omega_e= wpe - small noises. But, dt*wce > 1 is necessary for electron tracking.


### Execution Scripts ###

Linux: Compilation by mpif90, gfortran or PGI

 > mpich-4: ./configure --prefix=/opt/mpich-4 >&1 | tee conf.txt

 > fftw3: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/fftw3

mpif90 @mrg37-023A.f03 needs the parameter files param_A23A.h and rec_3d23A

 > $ mpif90 -mcmodel=medium -fast @mrg37-023A.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3

Execution by mpiexec (may need some tens of co-processors)

 > $ mpiexec -n 6 a.out &


### Simulation of Two Flux Bundles

One can enjoy simulations by changing system sizes and boundary conditions. For the present case, an equilibration of the pair of flux bundles is first tested in three dimensions. Fully kinetic ions and electrons are used, for example, in the rec_3d23A file. Then, let's start looking at a merging of two flux bundles. 
In-house graphic subroutines are incorporated in "@mrg37-023A.f03" in order to check the current run in the simulation. Figure 1 in the "EMfield.pdf" PDF plot shows the electric and magnetic fields in the YZ (left) and X (right) components at the early and final times. Two flux bundles at t=5000/wpe are seen touched and sqeezed at the y= Ly/2 plane. Reading papers of this implicit particle simulation code (Ref. 2-3) and applications to magnetospheric space plasmas (Ref. 1,4-5) are highly recommended.


### References: ###

1. M. Tanaka, Macro-particle simulations of collisionless magnetic reconnection, Phys.Plasmas, 2, 2920-2930 (1995).

2. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of three dimensions, J.Comput. Phys., 107, 124-145 (1993).

3. M. Tanaka, Macro-EM particle simulation method and a study of collisionless magnetic reconnection, Comput.Phys.Commun., 87, 117-138 (1995).

4. M. Tanaka, Asymmetry and thermal effects due to parallel motion of electrons in collisionless magnetic reconnection, Phys.Plasmas, 3, 4010-4017 (1996). 

5. H. Shimazu, M. Tanaka, and S. Machida, The behavior of heavy ions in collisionless parallel shocks generated by the solar wind and planetary plasma interactions, J.Geophys.Res., 101, 27565-27571 (1996).


