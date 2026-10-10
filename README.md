##  Water and Hydrate by TIP5P Molecular Dynamics Simulation ##

We have studied the microwave heating of water and ice by means of 
the TIP5P molecular dynamics simulation (with if_xyz2=.true.), 
and that of the ice and methane hydrate (with if_xyz1=.true.) by the same code.

This is the second generation code against the first one of 3-point SPC/E 
shake/rattle algorithm to study the microwave heating of water and ice. 
It was shown clearly in the first generation code that the ice below 273 K 
becomes a crystal and is not melted against the strong electric field 
of microwaves (JCP, 2007, Ref. 1).
The new results by the 5-body rotation coordinate system and TIP5P-Ewald model 
are shown below and the PDF in the arXiv Library: https://arxiv.org/abs/2311.01182 
(2023, 2024).

The five-point molecules for water are used known as the TIP5P-Ewald model. 
The five-water molecules are specified to calculate two hydrogens q_H= 0.241e 
and negative two hydrogens q_L= -0.241e with e the electron charge. 
The fifth oxygen atom of q_O= 0.0 is to correlate with adjacent molecules 
using the Lennard-Jones potential Psi(r)= eps_A/r^12 -eps_B /r^6 (Ref. 2).  

The fortran code with MPI is given in the file @p3mtip5p07a.f03 with additional files as param_tip5D07a.h, TIP507_config.start0, and initial coordinates 1cx666a.exyz and quaternions 1cx666a.q. Its description of the code is shown at README.md and also PDF documents of this repository. The latter documents are "Water_TIP5P_Simulation.pdf" for numerical coding, and for the physics simulation run by "Water_and_hydrate_molecules_by_TIP5P_code.pdf". 
The freezing ice state by microwaves, which is our theory discovery in JCP 2007 mentioned 
above, remains basically the same due to the structure of six-membered ice ! 

Methane hydrate is simulated by switching to if_xyz1=.true. of the TIP5P code, 
like with the SPC/E code (Ref. 3). We need the initial coordinate file mh3.exyz and quaternions mh3.q. 
It runs up to t= 1.7 x 10^(-8) s with E_x = 3. x 10^7 V/cm, which is terminated suddenly due to 
collapse of methane hydrate. Energies and scatter plots of molecules, the distribution of cosine's 
in the x-direction are shown in the PDF file "MDTip5pWater-8.pdf", the new version at https://arxiv.org/abs/2311.01182v8 (2024).

As natural resources, however, the methanes are 30 times environmentally more hazardous 
materials than carbon dioxide. The methane hydrate that will be mined and then burnt in air 
should be confined back to the deep interior of the earth if we need it as natural resources.


### Water Molecule Simulation by the TIP5P Model ###
   
   Motohiko Tanaka, Ph.D., Nagoya, Japan

Summary

a. Five sites are oxygen (O), hydrogens ($H_1, H_2$), and the virtual sites ($L_1, L_2$). They have charges of 0.0, 0.241e, -0.241e, respectively. The $L_1, L_2$ are  geometrically known from H-O-H and called the dummy sites.

b. Separate the translation $\bm{R}_{j}$, $\bm{V}_{j}$ and the rotation $\bm{r}_{i}$ for water where $j=1-N/5$ molecules and $i=1-N$ atoms, and $ (x_{i},y_{i},z_{i})$ means for the five sites. The separation is done at the starting step only; once determined at $t=0$, they are constant in time.

c. The half time step is first executed for a predictor, and the full step is made for a corrector.

d. Before the end of one step, the forces are calculated at $\bm{r}_{i}=\bm{R}_{j}+(A11*xr(i)+A21*yr(i)+A31*zr(i)), etc.$ with the three sites $i=1-3, 6-8, 11-13...$, and the L sites are calculated by algebraic operation. 

e. After correction of quaternions, go to the beginning of the cycle. The leap-frog method is used for the plasmas and waters.


* Each step consists of translation (1), rotation (2)-(4), and addition of field (5)-(8).

0. Read the quaternions from the file, read(30) e0,e1,e2,e3, j=1-N/5.

1. Sum the sites, $ d\bm{N_{j}}/dt=\sum_{k=1}^{5} \bm{F}_{i,k}/m_{j}, 
d\bm{R}_{j}/dt=\bm{V}_{j} $ with the translational motion. The (xr,yr,zr) and $ (Im_x,Im_y,Im_z) $ of molecules are the constants of time.

2. $ d\bm{L}_{j}/dt=(\sum_{k=1}^{5}
(y_{k}F_{i,k}^{z}-z_{k}F_{i.k}^{y},
 z_{k}F_{i,k}^{x}-x_{k}F_{i.k}^{z},
 x_{k}F_{i,k}^{y}-y_{k}F_{i.k}^{x}) $ for the
rotation matrix, summation over four sites.

3. $ \omega_{j,\alpha}=(A_{\alpha 1}L_{x}+A_{\alpha 2}L_{y}
+A_{\alpha 3}L_{z})/Im_{\alpha} $ for speices $A_{\alpha \beta}$
and moment inertia $\rm{Im}_{\alpha}$ for $\alpha=x,y,z$.

4. $ d\bm{q}_{j}/dt=Q(e_{i,0},e_{i,1},e_{i,2},e_{i,3})
(\omega_{j,x},\omega_{j,y},\omaga_{j,z}) $ for $\bm{q}$ of
Q and $\bm{\omega}$, by the Goldstein's book, and
get a new rotation matrix $$.

5. Get a new rotation matrix $ A(e_{0},e_{1},e_{2},e_{3}) $
of the time step.

6. \begin{equation*}
 \bm{r}_{i}=\bm{R}_{j}+
\begin{pmatrix}
A_{11} & A_{21} & A_{31} \\ 
A_{12} & A_{22} & A_{32} \\
A_{13} & A_{23} & A_{33} 
\end{pmatrix}
\begin{pmatrix}
xr_{i} \\ yr_{i} \\ zr_{i}
\end{pmatrix}
\end{equation*}
at three sites $\bm{r}_{i}$ plus the virtual sites (i=1-N), and $\bm{R}_{j}$ from \textit{Nr.1}.

7. Coulomb and LJ forces of the four sites and the oxygen site are calculated, which are  very heavily consumed of CPU. 

8. Corrections of quaternions at every 10 steps are executed, and go to the new time step.


### The Lennard-Jones Potential ###

With the Coulombic interactions, the Lennard-Jones 12-6 potential is adopted,
A=3.85x10^(-8) erg Ang^(12), B=4.36x10^(-11) erg Ang^(6) for TIP5P-Ewald sum's case, 
and A=4.17x10^(-8) erg Ang^(12), B=4.24x10^(-11) erg Ang^(6) for TIP4P.

Other parameters are: r(OH)=0.9572 Ang, Delta(HOH)=104.52 Ang. r(OM)=0.15 Ang is used 
for the TIP4P case, where the equipartition line of the virtual M sites is on the plane 
that equally separates the HOH angle. 

### In Order to Start a Run ###

To start a simulation of water cluster with the TIP5P code, the adjacent 4x4 hydrogen pairs of molecules 
are summed electrostatically, and also the oxygen pairs are coupled by TIP5P Lennard-Jones potentials (Ref. 5).

1. To get an initial state, we make the size of at least a 6x6x6 water cluster for numerical stability.
Short-range Coulombic and long-range forces are best separated for interactions, where the short-range 
forces for r>r_LJ ranges are made to be spatially dumped.

2. Around a given temperature, a dryrun is executed at least for 5 periods that is 10^(-9) seconds.
The long dryrun of a simulation is very important !!

3. Then, we may be ready to apply the electric field E_x= E_0*sin(\omega *time) in the x-direction to excite the 
electric dipole interactions of water. For the moment, we give the electric field 10 GHz where 
the electric field E_0 and electric dipole p_0 are of the order of 5x10^(-3) eV.
The E_0 value can be as small as 5x10^(-4) eV for the 6x6x6 water cluster.

### To Obtain the Initial Equilibrium for 298 K ###

To give uniform random noises to molecules, an equilibration takes extremely a long time of periods.
Thus, we try to use salt ions where 6x6x6 water clusters are fully strengthend with 64 Na(+) and 64 Cl(-) ions 
at t=1000, and waited by t=1700. Then, the salt ions are gradually removed by t=4700, and only water 
molecules are continued at least for 5 periods of t=50,000 before the real simulation run is again started from t=0.

The simulation run is shown in the next section in the PDF document with color pictures, 
"Molecular_Dynamics_Simulation_of_Water_and_Ice_by_TIP5P_Code.pdf".

It is noted that we impose the NVE ensemble simulation so that the box size is constant.
In the water heating of our daily experiment, however, the system is open and an energy goes out 
to the surrounding system. It is possible to have rapid heating of hot bubbling water.

### Figures of TIP5P-Model Simulations ###

Molecular dynamics simulations of water and ice, and/or methane hydrate, are shown by using the TIP5P-Ewald model, 
and the current file is "Water_and_hydrate_by_molecular_dynamics_TIP5P_code.pdf" in this repository (Ref. 5).

Figures in the color PDF file are water and ice and/or methane hydrate that include the energy of molecules, 
3-dimensional scatter plots, pair distribtion functions of oxygen and hydrogen atoms, and 
cosine distribution functions in the x-direction.

Generally speaking, water is heated by microwaves where its heating efficiency is highest at 0 Celsius.
Below 273 K, however, ice is frozen and not heated by the electric field.
Methane hydrate that is excited by microwaves application is heated and eventually collapsed. 

### Necessary Files of the Repository ###

Files for simulation runs: 

!*   1. @p3mtip5p07a.f03 : simulation code                       

!*   2. param_tip5p_D07a.h : parameter file 1, be constants      

!*   3. TIP507_config.start0 : parameter file 2, kstart=0 or 2   

!*      and/or TIP507_config.start1 : kstart=1 or 3              

!*   4. Initial files of water and ice: 1c666a.exyz and 1c666a.q

Post-processing files:                                

!*   * @lplotip507.f03 - the final history of energy             

!*   * @dipol_seqtip507.f03 - dipole Ex field                    

!*   * @iceplotip507.f03  - for 3D plot for x,y,z                

!*   * @wat_radtip507.f03 - radial distribution functions        
                                                         
References of Numerical Technique:

!* 1) M.Tanaka, J.Comput.Phys., vol. 79, 206 (1988).   

!* 2) M.Tanaka, J.Comput.Phys., vol.107, 124 (1993).   

!* 3) M.Tanaka, Comput.Phys.Comm., vol.87, 117 (1995). 

!*  4) M.Tanaka, Comput.Phys.Comm., vol.241, 56 (2019). 


### References ### 
 
1. M. Tanaka and M. Sato, Microwave heating of water, ice and saline solution: Molecular dynamics study, J.Chem.Phys., 126, 034509 (2007).

2. "Classical Mechanics", H. Goldstein, C. Poolee, J. Safko, 3rd Edition, Pearson Education Inc., England (2003); 古典力学，吉岡書店 (2006).

3. M. Tanaka, M. Sato, S. Nakatani, Microwave heating and collapse of methane hydrate by molecular dynamics simulations, https://arxiv.org/1909.01024, Cornell University, USA (2019).

4. M. Matsumoto, T. Yagasaki, and H. Tanaka,"GenIce: Hydrogen-Disordered Ice
Generator", J. Comput. Chem. 39, 61-64 (2017).

5. M. Tanaka, Molecular dynamics simulation of methane hydrate by means of the TIP5P-Ewald model,
https://arxiv.org/2311.01182, Cornell University Library, USA (2023, 2024).



