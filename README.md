Here is a repository of the necessary files to run the simulation of interacting 1-D particles, as implemented in my publications, using the Gillepsie Algorithm. For the deterministic mean-field equations using ODE time evolution, please see the accompanying repository: https://github.com/Blosberg/meanfield_void_numerics

A vido illustration of one such simulation is provided in the gif below; here we see red, triangular particles representing nucleosomes, and green squares, representing (smaller) transcription factors adsorbing on simulated promoter region with a (simplified) realistic energetic binding landscape. 

![](GA_GC_movie.gif) 

As such, nucleosomes are inhibited from binding in the NFR (100-200), and preferentially admitted to binding in the +1 regioni ( quasi-specifically bound near x~300), while TFs are preferentially bound in two arbitrarily chosen locations. Statistical positioning thereafter produces the familiar oscillatory pattern. A histogram of the time spent is shown below this region, while time progresses on a logarithmic scale.

The code-base itself can handle a much broader range of activity (e.g. active remodelling enzymes that grab and move particles around, arbitrary interaction potentials between neighbouring particles of different species, energetic landscape effects, etc.), and another version of this simulation was also developped for adsorption/desorption in 2D. These features are turned off by default, but can easily be reactivated.
 
If you would like to use this code for research purposes, please contact me. If you simply use it as is, I'd be content with a citation to the corresponding paper ( https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.088301 ), and perhaps an acknowledgement. If you would like me to adapt it to some specific use-case, or otherwise manage it for a project you are involved in, contact me, and we can discuss some arrangement for collaboration.

I've tried to add documentation here, but can't promise it will be perfectly "user-friendly" to a newcomer.

# Installation:

- To use this software, first clone this repository to your local machine, and then make sure that the GNU Scientific Library (GSL) is installed on your machine, and accessible on your "include" paths. You may need to adapt the Makefile that comes with this repo; I've tried to include paths to default directories for Homebrew. 

- Then, go to the repo directory on your machine and type `$ make `. If you're doing this on an OSX machine, you will likely get a compiler-dependent warning about dynamic vector sizing that you can safely ignore. On linux machines it should install smoothly.

- To run it, you should copy three files into your "Execution" directory. This is where all the data will come out, and you'll want to keep it separate from your scripts. The three files are the binary executable `GA_GC_N.x`, the input file `GA_GC_N.in`, and the seed for the random number generator `rngSEED.in`. Then navigate to your Execution directory and enter the following command, with 4 command line arguments:

` $ ./GA_GC_N.x  [NGtype]  [TASKID]  [muN] [eps]`

`[NGtype]` should have a value of either "HNG", "SNG", or "LNG", and describes the type of energetic interaction between neighbouring particles: "HNG" stands for "hard-core nucleosome gas" (i.e. mutual-exclusion), "SNG" stands for "soft-core nucleosome gas" (which implements the potential as described in https://www.pnas.org/content/110/14/5719), and "LNG" implements the "linear" potential, where interaction is directly proportional to proximity. 

`[TASKID]` should just be set to 1 in this case; this variable is used for batch submission on SGE clusters, so that different runs can implement different parameters.

`[muN]` provides the binding affinity of particles ("nucleosomes"), which can be realistically set to something around 12 for now.

`[eps]` describes the "stiffness" of particles, and is only relevant in the cases of "SNG" and "LNG". In the limit of infinite stiffness, the particles behave like "HNG" particles. See PRL for precise definition; a physical around where interesting dynamics take place would be in the range of `eps=20-25`, if `muN` is set to around 20. Again see PRL for ranges. 

# OUTPUT

Once the run is complete, a subdirectory using a name from the input file, and the input parameters will be created with various output results:
In particular, the "filling.." file shows density vs. time -- I recommend you plot time on a logarithmic scale.
Various other plots describe the profiles of gap-distribution between particles at various time points; the KL entropy; the 2-particle correlation between sizes of neighbouring gaps, and the overall frequency of occupation. For more details please get in touch.
