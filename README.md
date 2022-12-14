# GPU-accelerated classical Density Functional Theory (cDFT)
This repository archives related content for GPU-accelerated cDFT. The related content is still under development and for anyone interested in collaboration, please contact Jianzhong Wu (jwu@engr.ucr.edu) and Musen Zhou (mzhou035@ucr.edu)

©All rights reserved

## Try out code
The compiled file offers anyone who is interested in GPU-accelearated cDFT an oppoturnity to try it out. CUDA environment is required to excute both cDFT-MFA and cDFT-WDA. To run this code simply choose MFA or WDA with the example input provided, for example './MFA example.input'.

The first ouput is the signal from cg_descent package and 0 means convergence in cg_descent. the second output is the system pressure in the unit of bar and the following outputs are adsorption amount of gas molecules in the input in the unit of mol/L.


## Source Code
Source code of the compiled files can be found under the src/ folder.