from setuptools import setup, find_packages

setup(name='tetres',
      version='1.0.0',
      description='Time Tree Statistics package',
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
