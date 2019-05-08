from setuptools import setup
setup(
	name='pyDNM',
	version='0.1.0.0',
	url='https://github.com/dantaki/pydnm',
	author='Danny Antaki, Aojie Lian, James Guevara',
	author_email='dantaki@ucsd.edu',
	packages=['pyDNM'],
	package_dir={'pyDNM': 'pyDNM/'},
        package_data = {
            'pyDNM': ['*'],
        },
        scripts= ['pyDNM/pydnm'],
	requires=['numpy','sklearn', 'pandas']
)
