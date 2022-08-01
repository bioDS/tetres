from setuptools import setup, find_packages

setup(name='treeoclock',
      version='1.0.0',
      description='A package for discrete time trees',
      # url='http://github.com/storborg/funniest',  # todo
      author='Lars Berling',
      # author_email='flyingcircus@example.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['ete3',
                        'pandas',
                        'seaborn'
                        ],
      zip_safe=False)
