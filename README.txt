Authors: A. Taylor 		Universite catholique de Louvain,
         F. Glineur		Universite catholique de Louvain,
         J. Hendrickx	         Universite catholique de Louvain.

Date:   May 10, 2017

Version: May 10, 2017

----- Setup

In order to use the code, you should have YALMIP installed, along with some semidefinite solver (e.g. Mosek, Sedumi, SDPT3, ...).

Once YALMIP and the SDP solver installed (type 'yalmiptest' for checking the installation went fine); the toolbox can simply be installed
by executing the 'Install_PESTO' script (which only adds the required folders to your Matlab paths).
You can now execute the demo files for a step by step introduction to the toolbox.

----- Links

Link to YALMIP: http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download

Link to MOSEK: https://mosek.com/
Link to SeDuMi: http://sedumi.ie.lehigh.edu/
Link to SDPT3: http://www.math.cmu.edu/~reha/sdpt3.html

----- Approach

The toolbox implements the performance estimation approach as developped in the following articles:
 - "Smooth strongly convex interpolation and exact worst-case performance of first-order methods" (in Mathematical Programming).
 - "Exact Worst-case Performance of First-order Methods for Composite Convex Optimization" (Accepted in SIAM Journal on Optimization).
 
 Note that the approach of using semidefinite programming for obtaining worst-case guarantees was originally introduced in
 - Drori, Yoel, and Marc Teboulle. "Performance of first-order methods for smooth convex minimization: a novel approach." Mathematical Programming 145.1-2 (2014): 451-482.
 
----- Acknowledgments

The authors would like to thank Francois Gonze from UCLouvain and Yoel Drori from Google Inc. for their feedbacks on preliminary versions of the toolbox.


