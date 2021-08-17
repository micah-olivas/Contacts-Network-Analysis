from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
	name = 'contacts_network_analysis',
	version = '0.1.0',
	url = 'https://github.com/micah-olivas/Contacts-Network-Analysis.git',
	author = 'Micah Olivas',
	author_email = '',
	description = 'Contacts Network Analysis',
	packages = find_packages(),    
	install_requires = requirements,
)
