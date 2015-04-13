# DE_high-dimensional (Version 0.9.0)
=========================
Author: Ehsan Zahedinejad ezahedin@ucalgary.ca
=========================

Title:
=========================
Differential Evolution for high dimensional optimization problem.
See Publication here http://arxiv.org/abs/1501.04676
Program can be ran on single and mutiple processors.

Libraries:
============================================
MPI (Message Passing Interface), MKL (Math Kernal Library)
A single thread compilation of mkl leads to better efficiency. This is why I compiled using a single thread command.

Compilation
======================
Compilation command: ./make
However you should alway make sure that the program refers to the right location of libraries.
Edit makefile according to the location of you libraries.

