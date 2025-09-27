from setuptools import setup, find_packages

setup(
    name="polysurf",
    version="0.2.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "ase",
        "numpy",
        "scipy",
        "pandas",
        "matplotlib"
    ],
    python_requires=">=3.8",
    description="ASE-based builder to assemble slab + polymer systems for DFT workflows",
    url="https://github.com/ggbarbosa/polysurf",
    author="Gabriel Gouveia Barbosa",
    license="MIT",
)
