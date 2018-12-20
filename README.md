# Dynamic Occupancy Grid Maps (DOGMas)

Implimentation of the algorithm proposed in:

D. Nuss, S. Reuter, M. Thom, T. Yuan, G. Krehl, M. Maile, A. Gern, and K. Dietmayer. A
random finite set approach for dynamic occupancy grid maps with real-time application. arXiv,
abs/1605.02406, 2016.

The code demonstrates the particle filter processing of occupancy grids into DOGMas.

Run the 'run_all.sh' script from 'dogma/', which contains all the necessary commands in sequence. This sequence of commands generates a simple grid example from simulation and runs the particle filter velocity estimation.

The code is written in Python 2.7 with the following pip dependencies:

Cython==0.27.3

hickle==2.1.0

idna==2.6

ipykernel==4.7.0

ipython==5.5.0

ipython-genutils==0.2.0

ipywidgets==7.0.5

jupyter==1.0.0

jupyter-client==5.1.0

jupyter-console==5.2.0

jupyter-core==4.4.0

matplotlib==1.5.3

notebook==5.2.2

numpy==1.13.3

panda==0.3.1

pathlib2==2.3.0

pbr==3.1.1

pexpect==4.0.1

pickleshare==0.7.4

Pillow==4.3.0

scikit-image==0.13.1

scikit-learn==0.19.1

scipy==1.0.0

sklearn==0.0
