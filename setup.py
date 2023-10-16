#%% Setup the pip rfsed package
import setuptools

with open("./README.md", "r") as fh:
    long_description = fh.read()
REQUIRES = ['decorator', 'matplotlib>=2', 'numpy', 'scipy',
             'obspy>=1.0.3', 'tqdm', 'OWSLib',
            'cartopy', 'geographiclib', 'shapely', 'seaborn',
            'scipy', 'scikit-learn', 'h5','h5py', 'obspyh5', 
            'colorama==0.1.6', 'plotly==5.17.0', 'psutil==5.9.6']
CLASSIFIERS = [
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.9',
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
    author='Stephen Akinremi, Islam Fadel',
    author_email=['s.akinremi@utwente.nl', 'i.e.a.m.fadel@utwente.nl'],
    description="Package for receiver functions analysis and dealing with sediment effects in receiver functions",
    long_description=long_description,
    url="https://github.com/akinremisa/rfsed",
    packages=setuptools.find_packages(),
    install_requires=REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    classifiers=CLASSIFIERS,
    python_requires='>=3.6',
)
