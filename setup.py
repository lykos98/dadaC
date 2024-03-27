#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command, Extension
from setuptools.command.install import install

# Package meta-data.
NAME = "dadac"
DESCRIPTION = "C accelerated implementation of ADP"
URL = "https://github.com/lykos98/dadaC"
EMAIL = "francesco.tomba17@gmail.com"
AUTHOR = "Francesco Tomba"
REQUIRES_PYTHON = ">=3.6.0"
VERSION = "0.2.0"

# What packages are required for this module to be executed?
REQUIRED = [
    # 'requests', 'maya', 'records',
    "numpy",
    "scikit-learn",
    "wurlitzer",
]

# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}


EXT_DIR = os.path.join(os.path.dirname(__file__), "bin")



class RunMake(install):
    """Makefile on setuptools install."""

    def run(self):
        old_dir = os.getcwd()
        try:
            os.chdir(EXT_DIR)
            print("Building C library ...")

            os.system("make arm")
            os.system("make x86")
            #os.system("make lib")
            #self.spawn(['make'])

        finally:
            os.chdir(old_dir)
        install.run(self)


class RunMake_precompiled(build):
    """Makefile on setuptools install."""
    user_options = []


    def run(self):
        old_dir = os.getcwd()
        try:
            if(os.path.exists("bin/libdadac.so")):
                os.system("rm bin/libdadac.so")
            os.chdir(EXT_DIR)
            print("Building C library ...")
            os.system("make arm")
            os.system("make x86")
            #self.spawn(['make'])
        finally:
            os.chdir(old_dir)
        build.run(self)

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, "__version__.py")) as f:
        exec(f.read(), about)
else:
    about["__version__"] = VERSION


class UploadCommand(Command):
    """Support setup.py upload."""

    description = "Build and publish the package."
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print("\033[1m{0}\033[0m".format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status("Removing previous builds…")
            rmtree(os.path.join(here, "dist"))
        except OSError:
            pass

        self.status("Building Source and Wheel (universal) distribution…")
        os.system("{0} setup.py sdist bdist_wheel --universal".format(sys.executable))

        self.status("Uploading the package to PyPI via Twine…")
        os.system("twine upload dist/*")

        self.status("Pushing git tags…")
        os.system("git tag v{0}".format(about["__version__"]))
        os.system("git push --tags")

        sys.exit()

dadac_module = Extension(
    "dadac.core",
    sources = ["dadac/src/dadac.c", "dadac/src/kdtree.c", "dadac/src/kdtreeV2.c", "dadac/src/heap.c",  "dadac/src/vptree.c", "dadac/src/vptreeV2.c"],
    include_dirs=["dadac/include"],
    extra_compile_args=["-O3","-fopenmp"],
    extra_link_args=["-fopenmp", "-lm"]
)


# Where the magic happens:
setup(
    name=NAME,
    version=about["__version__"],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    #packages=find_packages(
    #    include=["dadac","dadac.src"],
    #    exclude=["tests", "*.tests", "*.tests.*", "tests.*","*src*","*include*","*src*"]
    #    ),
    packages=["dadac"],
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],
    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,

    license="MIT",
    package_data={"dadac": ["bin/*.so"]},
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    # $ setup.py publish support.
    #ext_modules=[dadac_module],
    cmdclass={
        'install': RunMake,
        'upload': UploadCommand,
        'build': RunMake_precompiled,
    },

)
