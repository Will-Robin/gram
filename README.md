# gram

Exploring the grammar of systems of chemical reactions by generating reaction
networks.

Calls to the chemoinformatics library are wrapped up into a single submodule so
that they can be easily replaced with a varying toolkits. Type annotations in
the code will help make sure that all of the components fit together correctly.

The current chemoinformatics toolkit in use is The RDKit. This package has
Python bindings, contains lots of functionality, and it widely employed in the
chemistry community. The network generation algorithms rely upon The RDKit's
reaction code. Performing chemical reactions in chemoinformatics packages is not
a trivial problem, so this code may behave imperfectly. However, it seems to
work well for simple organic reactions. Any of the core chemoinformatics code
can be switched out by changing the code in `gram/chemoinformatics/core` and
making sure that the input and output types are aligned in the
`gram/chemoinformatics` module.

For installation guide, see `installation.md`.

# Reaction Network Generation

A reaction rule is principally defined by reactive substructures and molecular
rearrangements. The reactive substructures component can be considered as a
'filter' for which compounds to choose to react with. This is an important
component of the reaction rule, which specifies which functional groups are
capable of reacting. The reaction transformation rule is responsible for making
sure that the transformations applied to molecules are done so with appropriate
accuracy (operating on the correct substructure of the molecule, and converting
it to the new substructure).

A reaction generation step runs as follows given (a) compound(s) and a reaction rule:

1. Check if the compounds have any substructures matching those contained in the
   reaction rule (find substructure matches in the reaction system).

   if not, pass

2. Apply the reaction rule to compounds with matching substructures.

