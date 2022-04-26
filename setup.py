from setuptools import find_packages, setup

setup(
    name="scope",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)
