# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
from pyepr import __version__, __copyright__
sys.path.insert(0, os.path.abspath('..'))

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

project = 'PyEPR'
copyright = __copyright__
author = 'Hugo Karas, Gunnar Jeschke, Stefan Stoll'
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.viewcode',
              'sphinx_design',
              'myst_parser',
              'sphinx.ext.intersphinx',
              'autoapi.extension',
              'sphinx_toolbox.collapse',
              'sphinx_toolbox.code',
              'sphinx_copybutton',
              'numpydoc',
              'sphinx_favicon',
              'sphinx_gallery.gen_gallery',
              'matplotlib.sphinxext.plot_directive',
              'sphinx.ext.imgmath'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for autoapi.extension -------------------------------------------
# https://sphinx-autoapi.readthedocs.io/en/latest/config.html
autoapi_dirs = ['../pyepr']
autodoc_typehints = "description"
autoapi_template_dir = "_templates/autoapi"

plot_include_source = True
plot_html_show_source_code = False
plot_formats = [('png', 100)]

autoapi_keep_files = True
autoapi_add_toctree_entry = False
autoapi_python_class_content= "both"
autoapi_python_use_implicit_namespaces = True
autoapi_own_page_level = 'class'

sphinx_gallery_conf = {
     'examples_dirs': 'examples',   # path to your example scripts
     'gallery_dirs': 'auto_examples',  # path to where to save gallery generated output
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_theme_options = {
    # "announcement": "Version 0.7 out now!",
    "light_logo": "pyepr_logo_light.svg",
    "dark_logo": "pyepr_logo_dark.svg",

}
favicons = [
    "favicon_16x16.png",
    "favicon_32x32.png",
    "icon.svg",
]