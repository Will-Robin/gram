# gram

A package for template-based reaction network generation in Python.

Calls to the chemoinformatics library are wrapped up into a single submodule so
that they can be easily replaced with a different chemoinformatics toolkit if
required. Type annotations in the code will help make sure that all of the
components fit together correctly.

For installation guide, see `installation.md`.

# Reaction Network Generation

A reaction rule is defined by reactive substructures and molecular
rearrangements applied to them. The reactive substructures component can be
considered as a 'filter' for which compounds to choose to react with, an
important component of the reaction rule, which specifies which functional
groups are capable of reacting. The reaction transformation rule is responsible
for making sure that the transformations applied to molecules are done so with
appropriate accuracy (operating on the correct substructure of the molecule,
and converting it to the new substructure).

A reaction generation step runs as follows given (a) compound(s) and a reaction rule:

1. Check if the compounds have any substructures matching those contained in the
   reaction rule (find substructure matches in the reaction system).

   if not, pass

2. Apply the reaction rule to compounds with matching substructures.

Examples of use are given in the `examples` folder.
