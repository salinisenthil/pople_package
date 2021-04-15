import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '21.4.1'
PACKAGE_NAME = 'pople'
AUTHOR = 'Salini Senthil, Raghunathan Ramakrishnan'
AUTHOR_EMAIL = 'salinis@tifrh.res.in'
URL = 'https://github.com/salinisenthil/pople_package'

LICENSE = 'MIT License'
DESCRIPTION = 'A toolkit for ab initio thermochemistry'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      include_package_data = True,
      install_requires=INSTALL_REQUIRES,
      #packages = ['.','templates'],
      packages=find_packages(),
      )

