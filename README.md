#  CS510 Robust Newton's Method Project
Latest Version: 0.1.0

---
This program takes an input polynomial and will iterate using any of four iterative fixed-point root-finding methods (Newton's method, Halley's method, Robust Newton's method, Hybrid Newton's method) on a region in the complex plane to visually display the convergence of points in the region to roots of the polynomial.  The actual roots' locations are also plotted.

The graphed polynomiograph is interactable!  One can zoom or pan as they please, and mouse over regions to see the index of the root to which they converge to.  The list of roots, as well as computation costs, are listed in the console.  
**NOTE: When using the .exe version, the console will close when the graph is closed.**

---
### Dependencies
* Python 3.7+
* Numpy
* Mathplotlib

---
### How to Use
Each version in the "Releases" tab above has a Python file `RobustVisualizer.py` and a Windows executable `RobustVisualizer.exe`. While either of these work, it is **highly** recommended that the Python file is run per the instructions below as it is much faster and more stable (and might not even be supported in the future).  Once could also run `main.py` straight from the source code, however this version is incredibly unstable and might not even work.

#### Instructions
1. Install Python 3 (verison 3.7 or newer) from [here](https://www.python.org/downloads/).
2. Run via command line:
	```
	python -m pip install --user numpy mathplotlib
	```
3. Within the same folder as `RobustVisualizer.py`, run via command line:
	```
	python main.py
	```
4. Follow the instructions given in the console window.

---
### Changelog
#### Version 0.1.0
* Initial release
* Supports 4 root-finding algorithms (Newton, Halley, Robust Newton, Hybrid Newton)
* Output is a 401x401 pixel graph (actual image size may vary) of the polynomiograph of the input polynomial
* All settings currently are hardcoded and can only be modified via the source file

---
### Credits
* SeaBass (smokedseabass@gmail.com) - Project Creator/Lead Programmer
* Bahman Kalantari (kalantar@cs.rutgers.edu) - Project Designer/Creator of "Robust Newton's Method"