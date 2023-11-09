import setuptools
from setuptools import setup

with open("./README.md", "r") as fh:
    long_description = fh.read()
REQUIRES = ['decorator', 'matplotlib', 'numpy', 'scipy',
             'obspy>=1.0.3', 'tqdm', 'OWSLib', 'pytest',
            'cartopy', 'geographiclib', 'shapely', 'seaborn',
            'scipy', 'scikit-learn', 'h5','h5py', 'obspyh5', 
            'colorama==0.1.6', 'plotly==5.17.0', 'psutil==5.9.6']
CLASSIFIERS = [
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Physics'
    ]
EXTRAS_REQUIRE = {
    'doc': ['sphinx', 'alabaster'], 
    'deconv_multitaper': ['mtspec'],
    'h5': ['obspyh5>=0.3'],
    'toeplitz': ['toeplitz'],
    'batch': ['tqdm']}
setuptools.setup(
    name="rfsed",
    version="0.0.1",
    author='Stephen Akinremi',
    author_email='s.akinremi@utwente.nl',
    description="Package for receiver functions analysis and dealing with sediment effects in receiver functions",
    long_description=long_description,
    url="https://github.com/akinremisa/rfsed",
    packages=setuptools.find_packages(),
    setup_requires=REQUIRES,
    install_requires=REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    classifiers=CLASSIFIERS,
    python_requires='>=3.6'
)