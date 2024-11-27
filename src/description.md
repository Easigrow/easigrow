# Easigrow

Summary
=======

Calculates fatigue crack growth rates and finds the optimum parameters for crack
growth models in order to match measurements for cracks grown under variable
amplitude spectrum loading.

Geometry
========

Geometry for a crack growth problem is specified via "sideways", "forward" and
"radius". The initial and final crack sizes are specified via “astart”,
“cstart”, ”aend” and “cend”. 

     a         crack length in thickness direction. 
     c         crack length in width direction. 
     forward   distance from crack origin to the nearest free edge in the
               thickness direction. 
     sideways  distance from crack origin to the nearest free edge in the width
               direction. 
     radius    hole or notch radius. 

Using the parameters above with a beta factor will define the geometry. 

     a) Surface crack in a plate         b) Corner crack in a plate

               |<-2c->|                        |<c>|                         
            ________________  ___               ________________  ___     
           |    \     / ^   |  ^               |   |   ^        |  ^      
           |     \___/  v a |  forward         |_ /    v a      |  forward
           |                |                  |                |         
           |________________| _v_              |________________| _v_     

                   |sideways|                  |<-- sideways -->|         


Usage
=====

Fatigue calculations and fatigue model optimisations are controlled from the
command line. Each command starts with “easigrow”. Think of each command as a
recipe, where you can add ingredients to: 

    1. Modify the fatigue sequence
    2. Count fatigue cycles
    3. Specify geometry
    4. Calculate fatigue crack growth
    5. Optimise model parameters
    6. Request for output

The commands to perform these tasks are listed below. By default Easigrow writes
fatigue crack growth calculation and model parameter optimisation results to the
console. This output can be redirected to a file using standard dos and linux
shell commands (e.g. “> out.txt”). In general, when a command requires a list
input, a comma-separated list with no spaces should be provided by the user.


Example
=======

easigrow --method Easigrow --seq_infile ft55.seq --scale 300 --seq_reorder
--cycle_method rainflow --beta seft-newman84  --dadn walker:default
--astart 10e-6 --cstart 20e-6 --aend 5e-3 --cend 10e-2 --output_vars blocks,a
--output_every 10  

This means: 

* Perform a fatigue crack growth prediction using the Easigrow method
  [--method Easigrow] 

* Read in the load sequence with the filename 'ft55.seq' [--seq_infile ft55.seq]
  and scale by a factor of 300 MPa [--scale 300] 

* Rotate the sequence so it starts and ends with the maximum peak [--seq_reorder],
  then extract fatigue cycles using the rainflow counting method
  [--cycle_method rainflow] 

* Calculate crack growth for a surface crack in a plate loaded in tension using
  the Newman and Raju beta solution [--beta seft-newman84]

* Start with a crack length in the a-direction = 0.01 mm and crack length in the
  c-direction = 0.02 mm. Terminate the calculation when a exceeds 5 mm or c
  exceeds 10 mm.

* Use the Walker da/dN vs. delta k crack growth model using the default
  parameters.  

* At the end of the calculation write out the crack length in the a-direction
  versus spectrum block, supplying output values every 10 blocks.
