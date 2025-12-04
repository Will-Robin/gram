"""
Chemoinformatic functionality module.

The code calling chemoinformatics libraries is in the `core` submodule so that
it can be interchanged without changing any of the other code.

To change the underlying chemoinformatics code, create a new submodule in the
`core` submodule. Name the functions inside the new code in a similar manner to
those in e.g. the rdkit submodule, and replace the rdkit import name in
`core/__init__.py`
"""
