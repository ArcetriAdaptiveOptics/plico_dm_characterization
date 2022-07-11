#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from shutil import rmtree

from setuptools import setup, Command

NAME = 'plico_dm_characterization'
DESCRIPTION = 'A suite of tools for calibration and characterization of deformable mirrors',
URL = 'https://github.com/lbusoni/plico_dm_characterization'
EMAIL = 'chiara.selmi@inaf.it',
AUTHOR = 'Chiara Selmi',
LICENSE = 'MIT'
KEYWORDS = 'Adaptive Optics, deformable mirror, plico, plico_dm, ALPAO, Boston MEMS, Physik Instrument',

here = os.path.abspath(os.path.dirname(__file__))
# Load the package's __version__.py module as a dictionary.
about = {}
with open(os.path.join(here, NAME, '__version__.py')) as f:
    exec(f.read(), about)


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


setup(name=NAME,
      description=DESCRIPTION,
      version=about['__version__'],
      classifiers=['Programming Language :: Python :: 3',
                   ],
      long_description=open('README.md').read(),
      url=URL,
      author_email=EMAIL,
      author=AUTHOR,
      license=LICENSE,
      keywords=KEYWORDS,
      packages=['plico_dm_characterization',
                'plico_dm_characterization.characterization',
                'plico_dm_characterization.configuration',
                'plico_dm_characterization.ground',
                'plico_dm_characterization.type',
                ],
      install_requires=["plico>=0.20",
                        "numpy",
                        "astropy",
                        "scipy",
                        ],
      include_package_data=True,
      test_suite='test',
      cmdclass={'upload': UploadCommand, },
      )
