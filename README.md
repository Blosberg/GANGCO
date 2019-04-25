Here is a repository of the necessary files to run the simulation of interacting 1-D particles, as implemented in my publications, using the Gillepsie Algorithm. 

A vido illustration of one such simulation is provided in the gif below; here we see red, triangular particles representing nucleosomes, and green squares, representing (smaller) transcription factors adsorbing on simulated promoter region with a (simplified) realistic energetic binding landscape. 

![](https://github.com/Blosberg/GA_GC/blob/master/GA_GC_movie_test2.gif) 

As such, nucleosomes are inhibited from binding in the NFR (100-200), and preferentially admitted to binding in the +1 regioni ( quasi-specifically bound near x~300), while TFs are preferentially bound in two arbitrarily chosen locations. Statistical positioning thereafter produces the familiar oscillatory pattern. A histogram of the time spent is shown below this region, while time progresses on a logarithmic scale.

The code-base itself can handle a much broader range of activity (e.g. active remodelling enzymes that grab and move particles around, arbitrary interaction potentials between neighbouring particles of different species, energetic landscape effects, etc.), and another version of this simulation was also developped for adsorption/desorption in 2D. 
 
If you would like to use this code for research purposes, please contact me, and  I would be happy to make it more "user friendly". At the moment, this repository has been left "as is" from production, and might not be easily incorporated into a newcomers workflow; I'd be  happy to add documentation and instructions if there were external interest to make it worthwhile.
