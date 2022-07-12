from setuptools import setup


setup(name='treeoclock',
      version='1.0.0',
      description='A package for discrete time trees',
      # url='http://github.com/storborg/funniest',  # todo
      author='Lars Berling',
      # author_email='flyingcircus@example.com',
      license='MIT',
      packages=['treeoclock'],
      install_requires=['ete3',
                        'pandas',
                        'seaborn'
                        ],
      zip_safe=False)
