import sys
try:
    import vtk
except ImportError:
    sys.exit("\"vtk\" is not properly installed. Please either install vtk and make sure \"import vtk\" wroks in python  or simply install vtkpytohn and run \"vtkpython setup.py install\". \nDownload vtk or vtkpython from: https://www.vtk.org/download/")

from setuptools import setup, find_packages

setup(name='molar',
      version='2.0.8',
      description='',
      url='https://github.com/alinar/Molar',
      author='Ali Narangifard',
      author_email='alinar@kth.se',
      license='GNU GENERAL PUBLIC LICENSE',
      packages=['molar'],
      install_requires = ['setuptools>=36.6.0','numpy'],
      zip_safe=False)
