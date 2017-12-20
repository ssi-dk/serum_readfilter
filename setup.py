from setuptools import setup, find_packages

setup(
    name='fastfilter',
    version='1.0',
    description='filters read against a small database of interest',
    url='https://github.com/ssi-dk/fastfilter',

    # Author details
    author='Kim Ng',
    author_email='kimn@ssi.dk',

    # Choose your license
    license='GPLv3',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3 :: Only',
    ],
    packages=find_packages(),
    scripts=['bin/fastfilter'],
    data_files=['config', ['fastfilter/config.yaml']],
    package_data={'fastfilter': ['test_data/*']},
    install_requires=[]
)
