Read me file for Fortran to HTML Converter

# Scripts for converting NAME III source code to HTML format.

This readme file gives brief instructions on running the scripts that
convert the Fortran 90 source code of NAME III into a HTML format for
viewing on a web browser.

+ Files in package:

 settings.txt       - contains information needed by the scripts
                      (the version number of NAME III and the input/output
                       directories must be set here)

 title_template.txt - a template .html file used for the title page

 keywords.txt       - a list of Fortran 90 keywords highlighted in the output

 info.sh            - preliminary script to create the cross-referencing
                      information (lists of functions, subroutines, modules
                      and types) needed by the main processing script

 fort2html.sh       - main processing script to convert Fortran 90 source code
                      into a web-based linked documentation system.

 readme.txt         - this readme file



+ Instructions:

 1) edit the three entries in the file "settings.txt" (the NAME III version
    number and the location of the input and output directories).
    IMPORTANT: although the order of the entries could be changed (if desired)
    the format of each line is FIXED - in particular, there must be NO white
    space around the "=".

 2) copy the .F90 files to be processed into the Input Directory.
    In principle, this could be the actual NAMEIII_Code directory of the
    model version being processed, but it's probably safer to create a
    duplicate of the files (I use a "ProcessingDirectory" adjacent to
    this Scripts directory for storing these temporary working files).

 3) run the preliminary processing script "info.sh" to create the
    cross-referencing information (files: functions.txt, subroutines.txt,
    modules.txt and types.txt).

 4) examine these four cross-referencing files to check that they are fairly
    sensible. The processing rules are rather basic and problems can sometimes
    occur here!
    IMPORTANT: the functions and subroutines associated with the Building
    module cause problems (because their names are often too general) and so
    it is safer to delete these entries in the cross-references.

 5) run the main processing script "fort2html.sh" to generate the HTML files.
 
 6) delete the files from the Input Directory after the script has completed,
    if required.
    

 + Known Problems:
 
  1) The functions and subroutines associated with the Building module
    can cause problems (because their names are often too general).
    This can be overcome by deleting their cross-references before running
    the main script. A longer term solution may be to standardise the
    Building module code to the style used in the rest of NAME III.
    
  2) The cross-referencing information (for subroutines, functions, types
    and modules) is created using rather simple search patterns. It is worth
    checking through the four cross-referencing files before running the main
    script in case of any problems.
    
    Please let me know of any further problems you may discover and I might
    be able to improve the pattern matching to address them.
