from distutils.core import setup

setup(
    name='destest',
    version='0.1',
    packages=['destest',],
    scripts=[],
    package_dir={'destest' : 'destest'},
    long_description=open('README.md').read(),
    )
