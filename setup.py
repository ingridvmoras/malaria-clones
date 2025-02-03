from setuptools import setup, find_packages

install_requires=[
    'numpy>=1.18.0',
    'pandas>=1.0.0',
    'matplotlib>=3.1.0',
    'seaborn>=0.10.0',
    'scipy>=1.4.0',
    'statannot>=0.2.3',
    'findpeaks>=1.0.5',
    'pytest>=5.0.0',
    'scikit-learn>=0.22.0',
    'matplotlib-venn'
]

setup(
    name='malaria-clones',
    version='1.0.0',
    description='A pipeline for peak detection and analysis in malaria clones data',
    author='Ingrid Vanessa Mora Sanchez',
    author_email='i.moras@uniandes.edu.com',
    url='https://github.com/ingridvmoras/malaria-clones',  
    packages=find_packages(),
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'peak-finder=peak_finder.main:main',  
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)