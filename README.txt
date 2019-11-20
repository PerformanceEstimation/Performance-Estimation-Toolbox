Authors: A. Taylor 		INRIA Paris,
         F. Glineur		Universite catholique de Louvain,
         J. Hendrickx		Universite catholique de Louvain.

Date:   Nov 2019

Version: Nov 2019

----- Introduction

(1) The document "PESTO_CDC2017_FINAL.pdf" contains the reference paper for the toolbox ("Performance Estimation Toolbox (PESTO): automated worst-case analysis of first-order optimization methods").
This paper contains a simplified and quick general introduction to the theory underlying the toolbox, and to its use.

If you used the toolbox, and found it useful for your research, you can mention the toolbox by citing this paper;
Taylor, Adrien B., Julien M. Hendrickx, and François Glineur. "Performance Estimation Toolbox (PESTO): automated worst-case analysis of first-order optimization methods." Proceedings of the 56th IEEE Conference on Decision and Control (CDC 2017). 2017.

(2) The document "UserGuide.pdf" contains more information and examples about the use of the toolbox.


----- Setup

In order to use the code, you should have YALMIP installed, along with some semidefinite solver (e.g. Mosek, Sedumi, SDPT3, ...).

Once YALMIP and the SDP solver installed (type 'yalmiptest' for checking the installation went fine); the toolbox can simply be installed by executing the 'Install_PESTO' script (which only adds the required folders to your Matlab paths).
You can now execute the demo files for a step by step introduction to the toolbox.

----- Links

Link to YALMIP: http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download

Link to MOSEK: https://mosek.com/
Link to SeDuMi: http://sedumi.ie.lehigh.edu/
Link to SDPT3: http://www.math.cmu.edu/~reha/sdpt3.html
