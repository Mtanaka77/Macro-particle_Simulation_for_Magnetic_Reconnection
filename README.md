## Macro-Particle Simulation for Magnetic Reconnection ## 

This page is discussed on the largescale electromagnetic particle simulation 
(J. Comp. Physics, Tanaka, 1993), which had a striking result of magnetic reconnection
in earth and astrophysical spaces (Phys. Plasmas, Tanaka, 1995, 1996). 
This was connected to heavy ions in collisionless parallel shocks (J.Geophys.Res., 
Shimazu, 1996)

### Magnetic Reconnection in Solar-Magnetospheric Couplings ###

Why is a large amount of the solar-eartth energy released in the distant magnetotail ?
This energy release is suddenly and typically observed as magnetic reconnection. 

There were many theories for the reconnection including from classical Dungey's theory to nuclear-fusion oriented anomalous resistivity. It is noted that Dr. Speicer paid attention as 'hypothesis' of inertia resistivity which was nearly forgotten as said to be unrealistic idea. 
Much later by a particle simulation as theory comptatinal tools, it was clearly shown that 'inertia of ions and electrons' has the key role of input and output flows for magnetic reconnection, resulting in large energy release of earth's magnetotail (Ref.1).

### Implicit and Explicit Particle-in-Cell Simulation Codes ###

An electromagnetic particle simulation code is utilized for solar and magnetospheric space physics (Ref. 1,4-5). 
The difference of the codes is that, for an explicit particle code, it is strictly bound by the Courant condition,  
Dx/Dt < c where Dx is the cell length, Dt is the time step, and c is the speed of light. 
On the other hand for the implicit particle code, it is free from this condition, that is Dx/Dt > c and is possible
to make research of solar physics environment. 

In the implicit case, both electric and magnetic fields are solved by the implicit condition 
where the low-frequency slightly backward time decentering technique is used. 
The backward decentering does not affect low frequency phenomena, \omega*Dt << 1 with
\omega = c/Dx (JCP, 1993, Ref. 2,3).
Magnetic reconnection and the solar wind-earth magnetic field coupling are quite suitable 
for applying this simulation code.

### Implicit Particle Simulation ###

One utilizes the time decentered scheme in \aimpl=0.6, while the time centered scheme in the explicit code (\aimpl=0.5) is used in other directory of molecular dynamics simulations. Four physical units are, i) time: 1/wpe (c/wpe: electron inertia length), ii) length: c/wpe, iii) mass: electron mass, and iv) charge: electron charge. The program is written in Fortran 2003 and is coded for parallelization by MPI ver.3.

The title, major references, and remarks of this simulation code are written in the top of the @mrg37_023A.f03 file.
Major subroutines are the following: /fulmov/ (particles are accumulated for their position and momentum), 
/emfild/ ans /cfpsol/ (electromagnetic fields are solved using large matrix equations), 
and /fulmov/(2) (particles are advanced), which are used in every time step. 

By the implicit scheme it is free from the Courant condition, that is, Dx(length)/Dt(time step) >< c, the speed of light. For the backward differential scheme in \aimpl > 0.5, a time step may be Dt~ 1.2/ \wpe in order to dump out plasma oscillations at plasma frequency \omega_e= \wpe - small noises. 

### Execution Scripts ###

Linux: Compilation by mpif90, gfortran or PGI

 > mpich-4: ./configure --prefix=/opt/mpich-4 >&1 | tee conf.txt

 > fftw3: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/fftw3

mpif90 @mrg37-023A.f03 needs the parameter files param_A23A.h and rec_3d23A

 > $ mpif90 -mcmodel=medium -fast @mrg37-023A.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3

Execution by mpiexec (may need some hundreds of processors)

 > $ mpiexec -n number_of_cpu a.out &


### Simulation of Two Flux Bundles

One can enjoy simulations by changing system sizes and boundary conditions. 
For the present case, an equilibration of the pair of flux bundles separating 
the poloidal magnetic field (the y-z component) is shown in three dimensions. 
Fully kinetic ions and electrons are used in the rec_3d23A file. 
Then, let's start looking at a merging of two flux bundles. 

In-house graphic subroutines are incorporated in "@mrg37-023A.f03" in order to check the current run in the simulation. 
Figure 1 in the first page of "EMfield.pdf" plot and Figure 2 in the following pages made at Nov. 2024 demonstrate 
the electric and magnetic fields of the macro-particle kinetic simulation.
The Y-Z (left) and X (right) components at the early and final times and  also 
precise plots of t= 0-4080/\wpe show merging of two flux bundles.
 Two flux bundles are seen touched and sqeezed at the y= Ly/2 plane. 
(However, it goes quite mixed and a simulation should stop at the final time.) 

Reading papers of this implicit particle simulation code (Ref. 2,3) and applications to magnetospheric space plasmas (Ref. 1,4,5) are highly recommended.


### References ###

1. M. Tanaka, Macro-particle simulations of collisionless magnetic reconnection, Phys.Plasmas, 2, 2920-2930 (1995).

2. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of three dimensions, J.Comput. Phys., 107, 124-145 (1993).

3. M. Tanaka, Macro-EM particle simulation method and a study of collisionless magnetic reconnection, Comput.Phys.Commun., 87, 117-138 (1995).

4. M. Tanaka, Asymmetry and thermal effects due to parallel motion of electrons in collisionless magnetic reconnection, Phys.Plasmas, 3, 4010-4017 (1996). 

5. H. Shimazu, M. Tanaka, and S. Machida, The behavior of heavy ions in collisionless parallel shocks generated by the solar wind and planetary plasma interactions, J.Geophys.Res., 101, 27565-27571 (1996).


