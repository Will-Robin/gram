from setuptools import setup, version

setup(
    name="gram",
    version="0.0.0",
    author="William E. Robinson",
    packages=["gram"],
    install_requires=["rdkit"],
    extras_require={"dev": ["pdoc", "black", "pylint"]},
)
