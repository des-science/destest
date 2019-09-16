from setuptools import setup

setup(
   name='destest',
   version='0.1',
   description='',
   license="MIT",
   author='Michael Troxel',
   author_email='michael.a.troxel@gmail.com',
   url="",
   packages=['destest'],
   install_requires=['fitsio', 'h5py', 'pickle', 'yaml', 'os', 'sys', 'time', 'cProfile', 'pstats', 'multiprocessing', 'matplotlib'],
)
