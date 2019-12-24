#  CS510 Robust Newton's Method Project
Latest Version: beta

---
This program takes an input polynomial and will iterate using any of four iterative fixed-point root-finding methods (Newton's method, Halley's method, Robust Newton's method, Hybrid Newton's method) on a region in the complex plane to visually display the convergence of points in the region to roots of the polynomial.

---
### Dependencies
* Python 3.7+
* Numpy
* Mathplotlib

---
### How to Use
Each version in the "Releases" tab above has a Python file `main.py` and a Windows executable `RobustVisualizer.exe`. While either of these work, it is *highly* recommended that the Python file is run per the instructions below as it is much faster and more stable (and might not even be supported in the future).  Once could also run the `main.py` straight from source, however this version is unstable and might not even work.

#### Instructions
1. Install Python 3 (verison 3.7 or newer) from [here](https://www.python.org/downloads/).
2. Run via command line:
    ```
    python -m pip install --user numpy mathplotlib
    ```
3. Within the same folder as `main.py`, run via command line:
    ```
    python main.py
    ```

---
### Credits
* SeaBass (smokedseabass@gmail.com) - Project Creator/Lead Programmer
* Bahman Kalantari (kalantar@cs.rutgers.edu) - Project Designer/Creator of "Robust Newton's Method"