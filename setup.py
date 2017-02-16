from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='hades',
      version='0.1',
      description='fluid mechanics, aero and post processing tools',
      url='',
      author='ISAE/DAEP',
      author_email='',
      license='',
      packages=['hades'],
      zip_safe=False)
