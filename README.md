bloch-simulator-python
======================

The original bloch equation simulator was a Matlab mex file created by Brian Hargreaves at Stanford University. This is a modification to run it as a Python C extension
We used the simulator in a graduate MRI class taught by Mikki Lustig; Lustig wrote several helper modules in matlab, which I've also converted to Python.
This module current uses python3. If there is demand for it, I can write a version for python2 as well.

Dependencies
======================
python 3.4
numpy 1.8.1
scipy 0.14.0
matplotlib 1.3.1

These are the version of the libraries I have tested the sim on. I make no guarentee about other versions, but I am not using very complicated calls, so there should be some flexibility.

Installation
======================
Simply run "python setup.py install" to install the simulator. Then "from bloch import bloch" for the primary bloch simulator function.

License
======================
This library is distributed under the same terms as Brian's original bloch simulator

Thank you to Brian for the original bloch simulator, Mikki Lustig for the code for the helper modules, and NPann for assisting me with a critical bug fix.
