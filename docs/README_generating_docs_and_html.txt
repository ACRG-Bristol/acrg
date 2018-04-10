
+++++++
Summary
+++++++

To update the docs and generate new HTML, run the following commands:
>> sphinx-apidoc -f -o api $ACRG_PATH
>> make html

To look at the generated HTML open _build/html/index.html with any internet browser e.g.
>> firefox _build/html/index.html

Sphinx can be used directly with a documentation hosting service such as ReadTheDocs if the repository
is public (and so can be seen by ReadTheDocs to create the documentation). In the meantime, we can keep 
this maintained as a nice format for collecting together all of our documentation into a centralised set
of web pages.

--------
Creation
--------

This documentation folder has been generated using Sphinx using the command:
>> sphinx-quickstart

Questions are posed on the command line about the setup including:
 - Root path for the documentation: docs/
 - Project name: acrg
 - Project author: mrghg
 - Please indicate if you want to use one of the following Sphinx extensions:
   - autodoc: automatically insert docstrings from modules (y/n) [n]: y
   - todo: write "todo" entries that can be shown or hidden on build (y/n) [n]: y
   - coverage: checks for documentation coverage (y/n) [n]: y
   - imgmath: include math, rendered as PNG or SVG images (y/n) [n]: y
 - A Makefile and a Windows command file can be generated for you so that you only have to run e.g. 
 `make html' instead of invoking sphinx-build directly.
   - Create Makefile? (y/n) [y]: 
   - Create Windows command file? (y/n) [y]: 

---------------
rst and autodoc
---------------

Sphinx uses reStructured text files (.rst) as the basis of its documentation model.

See docutils.sourceforge.net/docs/ref/rst/restructuredtext.html for more details of this format.

These can be created and edited directly but can also be generated from the docstrings within our module
files (contained within ./api folder) using autodoc. * Unless using some external hosting service, such
as ReadTheDocs, this will not automatically update when we make changes to the scripts *.
As described in the Summary section, these files are created/updated using the sphinx-apidoc command.

For more details, see: http://www.sphinx-doc.org/en/master/man/sphinx-apidoc.html

To see how this would look, as described in the summary section, the make html command can be used to
show what the html would look like.

