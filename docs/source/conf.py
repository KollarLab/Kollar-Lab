# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os, sys
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))
sys.path.append('C:/Program Files/IVI Foundation/IVI/Bin')
sys.path.append('C:/Program Files/Acqiris/MD3/bin')


print(sys.path)

project = 'Kollar-Lab Repository'
copyright = '2024, Kollar Lab'
author = 'Kollar Lab'
release = '2024'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              ]


autosummary_generate = True
#add_module_names = False

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

html_theme_options = {
    "show_nav_level": 0,
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink"],
    "navbar_align": "right",
}

html_sidebars = {
    '**': [
        'localtoc.html',
        'globaltoc.html'
    ]
}