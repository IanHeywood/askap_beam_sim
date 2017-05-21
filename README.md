# askap_beam_sim

Implementation of parametrised fits to the holography beams as 2D Gaussians inside a MeqTrees simulation.

askap-sim.py is a stripped down version of turbo-sim.py from MeqTrees' [Siamese framework](https://github.com/ska-sa/meqtrees-cattery/tree/master/Cattery/Siamese) which invokes askap_beams.py to simulate the beams. 

The 36 beams are characterised in askap12_beams.p which a single dictionary, with keys of the form ANT:BEAM:POL and value tuples of (x0,y0,bmaj,bmin,pa,A).

Since antenna 2 was used as a reference antenna for the holography that measured this set of beams it has been assigned values in that tuple that are the medians of the per-beam values for the other 11 antennas.
