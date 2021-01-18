from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='aerokit',
      version='0.1',
      description='fluid mechanics, aero tools',
      url='',
      author='ISAE/DAEP',
      author_email='',
      license='',
      packages=['aerokit'],
      zip_safe=False)
