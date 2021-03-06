from setuptools import setup, find_packages

setup(
    name='serum_readfilter',
    version='1.0',
    description='filters read against a small database of interest',
    url='https://github.com/ssi-dk/serum_readfilter',

    # Author details
    author='Kim Ng',
    author_email='kimn@ssi.dk',

    # Choose your license
    license='GPLv3',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        'Natural Language :: English',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'serum_readfilter = serum_readfilter.serum_readfilter:main',
        ]
    },
    include_package_data=True,
    data_files=[('config', ['serum_readfilter/config.yaml'])],
    package_data={'serum_readfilter': ['tests/*']},

    test_suite='nose.collector',
    tests_require=['SerumReadFilter']
)
