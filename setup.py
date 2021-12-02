from setuptools import setup

setup(name='treeoclock',
      version='0.0.1',
      description='A package for time trees',
      # url='http://github.com/storborg/funniest',
      author='Lars Berling, Lena Collienne, Jordan Kettles',
      # author_email='flyingcircus@example.com',
      license='MIT',
      packages=['treeoclock'],
      install_requires=['ete3',
                        'pandas',
                        'seaborn'
                        ],
      zip_safe=False)