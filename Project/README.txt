WARNING:
This code was written by a lazy undergrad.
Do not use this code for anything that even somewhat matters.
This code is slow and probably buggy (though i'm not aware of any bugs).



To use this code, make a *.py file that imports diffusion_solver.
Looking at the sample code in demo.py is probably more helpful than reading this file.


The diffusion_solver class requires a 2D array of PhysicsVals with the following properties:

D: the diffusion coefficient
neutrons travel a long distance in a material with a large D and a short distance in a material with a small D
assigning D=0 means that the neutrons will simply sit still and pile up until equilibrium (flux = S/Sigma_a)
assigning D<0 means that neutrons travel against the concentration gradient which is impossible and
does not converge and will cause an error

Sigma_a: the absorption cross section
materials with large Sigma_a will absorb more of the neutrons in the material than those with small Sigma_a
assigning Sigma_a = 0 in some places is perfectly fine
assigning Sigma_a = 0 everwhere with reflecting boundaries on all sides means the neutron population will
increase without bound
Sigma_a < 0 implies that the material is a multiplying medium, that is, a neutron in the material will produce
additional neutrons (e.g, through fission)
This code is not capable of doing criticality calculations, so if you assign a very negative Sigma_a, the neutron
population will multiply out to infinity and cause an error
If Sigma_a is negative but the system is subcritical,
this code will solve for the subcritical equilibrium neutron concentration, but this hasn't been tested

S: neutron source
Cells with nonzero S will be an independent source of neutrons
Making S negative is nonphysical, but won't cause an error, in case you want to use this code for something that it
makes sense to have a negative quantity of. This will likely give bad results. See warning at top.


The diffusion_solver class requires arrays giving the widths of each column and heights of each row in the mesh.
These should all be greater than zero.

You also need to supply boundary conditions.
Assigning a number to a boundary condition will fix the flux at the boundary to the number. (Dirichlet boundary condition)
Assiging 'ref' makes the boundary a reflecting boundary.
Other boundary conditions are not supported.

Once all of these things are assigned to the diffusion_solver, simply call its solve() method.
For a 128 by 128 cell mesh, this takes 10-15 seconds on a 2012 laptop
(longer with especially unreasonable physics values) and scales as O(N^3), so don't try to
solve a very large system using this code. (See also warning at top)
Memory use scales as O(M*N) or O(N^2) for a square mesh

To execute the code, do to the *.py file you made whatever it is you normally do to run python programs, e.g. write

python demo.py

in the terminal


tests.py contains only four tests:

compareAnalytical1DFixed compares the result of the diffusion solver to an analytical solution
of the diffusion equation in one dimension for vacuum boundaries and constant physics values

compareAnalytical1DCos is the same, but uses a source distribution that varies as cos(x)

compareAnalytical2DFixed compares the result with a 2D fixed source solution for vacuum boundaries on all four sides

verifyReflection ensures that reflecting boundaries are implemented correctly by solving the same system twice;
once using reflecting boundary conditions, and once by actually reflecting the mesh and solving a mesh
four times as large with fixed flux boundary conditions on all sides

all use randomized values and were used to generate cosh.png, cos.png, cosh2.png, and rorschach.png, respectively