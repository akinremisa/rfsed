import os
import sys

# Add your package's root directory to sys.path here.
# For example, if the package is one level up:
sys.path.insert(0, os.path.abspath('../src'))

# Import the package to get the version
import rfsed

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'rfsed'
copyright = '2024, Stephen Akinremi'
author = 'Stephen Akinremi'
release = rfsed.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']

autodoc_member_order = 'bysource'

html_logo = './logo/rfsed_logo_name.png'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_theme = 'alabaster'
# html_static_path = ['_static']
