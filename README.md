# bsc-code
The code for the Bachelor of Science project, performed at the University of Cologne under the supervision of Achim Rosch.

Abstract:
In this theoretical thesis we make a study of the Coulomb disorder in the bulk of a
topological insulator (TI) nanowire, for example, made of (Bi_{1−x} Sb_x )_2 Te_3 . As many
topological insulator materials are doped semiconductors, randomly placed charged
impurities in their bulk create disordered electrostatic potential, whose fluctuations
may result in formation of bulk charge puddles, possibly having an effect on trans-
port properties of these materials.
In the first part of this work (sections 3 to 5), we investigate the surface screen-
ing of a TI nanowire, introducing a numerical algorithm which allows for making
computations of surface potential and charge density for various conduction regimes
depending on the surface chemical potential.
In the second and final part of this work (section 6), we simulate Coulomb disorder
in the nanowire bulk including the effect of surface screening, in order to analyse
the possibility of bulk puddle formation dependent on the strength of the surface
chemical potential and the width of the bulk band gap. For that, we study two key
quantities: the neutral dopant density and the screening length scale close to the
nanowire surface.


The code is organised in a number of scripts, including different useful functions. The heart of the code is the Efros-Schklovskii algorithm, written in the file efcAlgo.jl.
