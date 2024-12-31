from setuptools import setup, find_packages

setup(
    name='malaria_clones',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'seaborn',
        'scipy',
        'statannot',
        'findpeaks',
        'pytest',
        'scikit-learn'
         
    ],
    entry_points={
        'console_scripts': [
            'malaria_clones=malaria_clones.main:main',
        ],
    },
    author='Ingrid Vanessa Mora Sanchez',
    author_email='i.moras@uniandes.edu.co',
    description='A project to identify peaks of malaria parasite density and count first appearences of clones',
    url='https://github.com/ingridvmoras/malaria_clones',
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
)
