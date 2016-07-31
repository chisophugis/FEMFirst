This is the beginning of a finite element solver.

The problem that it solves is a discrete roller-spring system as
described in Bathe's lecture
https://www.youtube.com/watch?v=oNqSzzycRhw
(you could of course use it to solve analogous problems like resistive
networks).

The code is in `main.cpp`.
The file `input.txt` shows an example of the inputs.


I spent a lot of time studying MFEM (http://mfem.org/) for guidance.
MFEM is really well written and tremendously helped me understand FEM.
The OCW course that the above Bathe lecture belongs to is quite good
and also helped a lot.


The routine `cg` is a conjugate gradient solver.
One of the main reasources for writing this was Shewchuk's "An
Introduction to the Conjugate Gradient Method Without the Agonizing
Pain".
The code in `cg` is not just blindly copied out of a book/paper.
My derivations are in `cgderivation/` (of course done with a bunch of
staring at other resources).


The routine `eliminateBoundaryConditions` is actually one of the parts
that I found most difficult to write, as I couldn't find any resource
that actually described the intuition for how to do what it does
(hence the extensive comments).
Most of my understanding comes from staring at MFEM.

