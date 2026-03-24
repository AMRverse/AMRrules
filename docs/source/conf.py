# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AMRrules'
copyright = '2026, Jane Hawkey, Kathryn Holt, Natacha Couto, ESGEM-AMR Working Group'
author = 'Jane Hawkey'
release = '1.0.0'
version = '1.0.0'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_design",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_theme_options = {
    "show_toc_level": 3,
    "home_page_in_toc": True,
    "navigation_depth": 4,
    "repository_url": "https://github.com/AMRverse/AMRrules/",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": False,
    "path_to_docs": "docs", # Change this if your folder is named 'doc' or similar
    "repository_branch": "main",

}
#html_static_path = ['_static']
html_logo = 'https://github.com/AMRverse/AMRrules/blob/genome_summary_report_dev/docs/source/images/AMRrules_logo.png?raw=true'
#html_favicon = 'https://github.com/klebgenomics/Kaptive-Web/blob/master/static/images/favicon.png?raw=true'
