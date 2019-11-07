from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='atmopy',
      version='2.0.33',
      description='Package to plot and process atmospheric data',
      url='https://dargueso@bitbucket.org/dargueso/atmopy.git',
      author='Daniel Argueso',
      author_email='d.argueso@uib.es',
      license='MIT',
      packages=['atmopy'],
      install_requires=[
          'numpy','netcdf4','wrf-python','datetime','pandas','matplotlib',
      ],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['bin/get_wrfvar.py']
)
