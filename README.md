# gram

Exploring the grammar of systems of chemical reactions.

Adapted from NorthNet code; stripped down and rearranged.

The chemoinformatics libraries have been wrapped up into a single submodule so
that they can be easily replaced. Type annotations in the code will help make
sure that all of the components fit together correctly.

# Installation

1. Clone the repository and navigate into it:

```shell
git clone https://github.com/Will-Robin/gram.git
cd gram
```

2. Create a virtual environment via your preferred method and activate it. For
   example, using `conda`:

```shell
conda create -n gram-env
conda activate gram-env
```

3. Install gram in the environment:

```shell
pip install -e .
```

4. Test whether you can call `gram` in Python scripts (the following code should
   not give an error:

```shell
echo from gram import Classes > test.py
python test.py
```

(then remove `test.py` with `rm test.py`)
