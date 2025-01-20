from setuptools import setup, find_packages

setup(
    name='malaria-clones',
    version='1.0.0',
    description='A pipeline for peak detection and analysis in malaria clones data',
    author='Ingrid Vanessa Mora Sanchez',
    author_email='i.moras@uniandese.edu.com',
    url='https://github.com/ingridvmoras/malaria-clones',  # Replace with your project's URL
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'seaborn',
        'matplotlib',
        
    ],
    entry_points={
        'console_scripts': [
            'run-pipeline=main:main',  
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)