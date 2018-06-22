# LIGGGHTS-NANO

This software package is based on LIGGGHTS-PUBLIC v3.8 developed by DCS computing (www.cfdem.com). 

This version includes contact models for aggregated nanoparticles smaller than 25 nm. At this scale, the dominating capillary and solvation forces must be considered. Details can be found in: 
https://link.springer.com/article/10.1007/s10035-018-0799-9

another article (considering the bond model) is currently under review and will be added later.

All new files can be found in src/Nano:
- Files necessary to reimplement the atom style molecular (atom_vec_molecular.h/.cpp)
- Files for contact models:
    + cohesion_model_bond.h
    + cohesion_model_bond_stiffness.h
    + normal_model_hertz_stiffness_nano.h
    + rolling_model_nano.h
    + tangential_model_nano.h
- Files for tabulated potentials (pair_table.h/.cpp)
- multi_node_mesh_parallel_I.h, which was adjusted to allow walls with dimensions smaller than 1e-4 m at unit style si

## Some information on the bond model
The applied bond model is based on Potyondy and Cundall: A bonded-particle model for rock (2004) https://www.sciencedirect.com/science/article/pii/S1365160904002874

It portation of LIGGGHTS-WITH-BONDS (https://github.com/richti83/LIGGGHTS-WITH-BONDS), which uses the LAMMPS bond interface.
This version uses a cohesion to implement bonds. 

It is further optimized to our purpose. It creates a bond at timestep 0 (initial timestep) and has no option to recreate broken bonds or to have another trigger for bond creation. Since it was not in our scope, we did not try to write a fix that controls this behavior. 
It further has the option that bonds are only created between particles of the same molecule (that's why the molecular atom style is reimplemented). 
This is controlled by the *pair_style* command. Just add *per_molecule on* or *per_molecule off* at the end of the line. Default is *on*.

There are two versions implemented: *bond* and *bond/stiffness*.
*bond/stiffness* is the version of LIGGGHTS-WITH-BONDS. We preferred, however, to calculate the stiffness in runtime, to account for flexible bond lengths. 
Therefore: If *bond* is used: you can parse the option *full_length on/off* at the end of the *pair_style* line. If it is set to *off*, the distance from particle surface to particle surface is used for the bond length (distance of centers of mass, otherwise).

## Calculation of the tabulated potential
The capillary and solvation force is highly computationally demanding. Therefore, we have implemented them via tabulated potentials. 
We provide a tool that generates the tabulated potentials and also calculates the contact stiffness *kn*/*kt*. You can find it in the folder *writePotentialFilesNano/generatePotentialFiles.py*. It requires *scipy* for integration (the potential tables provide the energy additional to the force)

within *generatePotentialFiles.py* you can define the diameters you need and some other values. It automatically generates all necessay files (*include.potential_table* and *include.potential_params*). Please note, that the bond parameters are not included in these files. You need to parse them independently.

## Selection of the bond models
The following lines enable the selection of all models:  

pair_style      hybrid/overlay gran model hertz/stiffness/nano tangential nano cohesion bond rolling_friction nano per_molecule on table linear 450  
pair_coeff        * * gran  
pair_coeff        1 1 table include.potential_table CAP_SOLV_3.0_3.0  
pair_coeff        1 2 table include.potential_table CAP_SOLV_3.0_4.0  
pair_coeff        1 3 table include.potential_table CAP_SOLV_3.0_5.0  
pair_coeff        1 4 table include.potential_table CAP_SOLV_3.0_6.0  
pair_coeff        1 5 table include.potential_table CAP_SOLV_3.0_7.0  
pair_coeff        1 6 table include.potential_table CAP_SOLV_3.0_8.0  
pair_coeff        1 7 table include.potential_table CAP_SOLV_3.0_9.0  
pair_coeff        1 8 table include.potential_table CAP_SOLV_3.0_10.0  
pair_coeff        1 9 table include.potential_table CAP_SOLV_3.0_11.0  
pair_coeff        1 10 table include.potential_table CAP_SOLV_3.0_12.0  
pair_coeff        1 11 table include.potential_table CAP_SOLV_3.0_13.0  
pair_coeff        1 12 table include.potential_table CAP_SOLV_3.0_14.0  
pair_coeff        1 13 table include.potential_table CAP_SOLV_3.0_15.0  
pair_coeff        1 14 table include.potential_table CAP_SOLV_3.0_16.0  
.  
.  
.  

pair_style defines the contact models as usualy. The solvation/capillary force model is superimposed with the gran model (hertz/stiffness/nano model). Therefore, hybrid/overlay is applied to combine these two. 
**pair_coeff      ** gran** enables the gran model (hertz/stiffness/nano) for all particle particle interactions. The following lines **pair_coeff** define the tabulated potentials per atom_type combination. e.g.:
type 1 - type 1 is taken from include.potential_table and applies the table CAP_SOLV_3.0_3.0 for the interactions. To see the structures of the tabulated potentials, see the example folder.

## Here is a list of required material parameters for the models:

**fix     p1 all property/global kn peratomtypepair**  
This is the normal contact stiffness in N/m^2

**fix     p2 all property/global kt peratomtypepair**  
This is the tangential contact stiffness in N/m^2 (usually equal to kn)

**fix     p3 all property/global coefficientRollingFriction peratomtypepair**  
The rolling friction coefficient N/m^2. The friction is dependent on contact area, therefore, the rolling friction coefficient is not implemented in the classical sense. This is obviously confusing and might be corrected later

**fix     p4 all property/global coefficientFriction peratomtypepair**  
The tangential friction coefficient N/m^2. The friction is dependent on contact area, therefore, the rolling friction coefficient is not implemented in the classical sense. This is obviously confusing and might be corrected later.

**fix     p5 all property/global gamman peratomtypepair**  
The damping coefficient for normal contact in Ns/m^3. This will be multiplied with the sqauare of the effective particle radius (during runtime)

**fix     p6 all property/global gammat peratomtypepair**  
The damping coefficent for tangential contact. This is set to zero for the current model. It needs to be parsed, however.

### Bond Model **bond**

**fix     p7 all property/global bondYoungs peratomtypepair**  
Youngs modulus of the bond (sinter bridge) in Pa

**fix     p8 all property/global bondPoisson peratomtypepair**  
Poission's ratio of the bond (sinter bridge)

**fix     p9 all property/global bondTensileStrength peratomtypepair**  
Tensile strength of a bond in Pa

**fix     p10 all property/global bondShearStrength peratomtypepair**  
Shear strength of a bond in Pa

**fix     p11 all property/global bondRadius peratomtypepair**  
The radius of the bond. It is typically lambda * min(R1,P2), with lambda being a constant value. In the example, 0.48 is used

**fix     p12 all property/global bondMaxSeparationDistance peratomtypepair**  
Threshold distance between the particles in the initial timestep

### Bond Model **bond/stiffness**

**fix     p7 all property/global bondSn peratomtypepair**  
Bond normal stiffness

**fix     p8 all property/global bondSt peratomtypepair**  
Bond tangential stiffness

**fix     p9 all property/global bondTensileStrength peratomtypepair**  
Tensile strength of a bond in Pa

**fix     p10 all property/global bondShearStrength peratomtypepair**  
Shear strength of a bond in Pa

**fix     p11 all property/global bondRf peratomtypepair**   
This factor is multiplied with min(R1,R2) to determine the bond radius

**fix     p12 all property/global bondMaxSeparationDistance peratomtypepair**  
Threshold distance between the particles in the initial timestep


If there are any questions regarding the models or the implementation, please do not hesitate to contact Valentin Baric (v.baric@iwt.uni-bremen.de) or Prof. Lutz Mädler (lmaedler@iwt.uni-bremen.de). Please refer to the papers given at the top of the document.
