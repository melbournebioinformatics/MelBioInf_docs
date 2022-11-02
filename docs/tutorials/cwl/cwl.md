
<img src='https://www.melbournebioinformatics.org.au/wp-content/themes/_vlsci/images/logos.png' style="width:350px; padding-right:50px"><img src='https://d2glwx35mhbfwf.cloudfront.net/v13.3.2/logo-with-padding.svg' style="width:150px">

# Introduction to CWL

>**Important**<br>
>This material is no longer maintained and may be out of date.<br>
>Please go to this link below for a current [Introduction to CWL](../cwl_2022/cwl_2022.md)

## Summary
Common Workflow Language (or CWL), is a growing language for defining workflows in a cross-platform and cross-domain manner. In biology in particular, we need workflows to automate complex analyses such as DNA variant calling, RNA sequencing, and genome assembly. CWL provides a simple and well-defined format for automating these analysis by specifying their stages and connections using readable CWL documents.
 
CWL makes use of a number of existing standards, including support for cluster computing using SLURM or PBS, containerisation using Docker, and deployment using common packaging formats. In addition, the CWL ecosystem has grown to include workflow visualisation tools, graphical workflow editors, libraries for interacting with CWL programatically and tools that convert to and from CWL and other workflow formats.

## Outcomes
At the end of the course, you will be able to:

* Find and use CWL tool definitions online
* Use the Rabix Composer, a graphical editor for CWL
* Understand how to write CWL tool definitions for command line tools
* Use Docker with CWL to provide software dependencies and ensure reproducibility
* Join CWL tools into a workflow
* Read and write CWL files written in YAML
* Understand advanced CWL features like secondary files, parameter references and subworkflows
* Run CWL workflows on local and HPC systems 
 
## Requirements
### General
* This workshop is aimed at anyone with basic Unix command-line experience
* Attendees are required to bring their own laptop computers

### Software

All of this software is free, and should run on any operating system (Mac, Windows, or Linux):

* Rabix Composer: 
    * A GUI for writing CWL
    * <https://github.com/rabix/composer/releases>
    * Download the .dmg (Mac), .exe (Windows) or .AppImage (Linux)
* Docker
    * A system for managing containers
    * <https://store.docker.com/search?type=edition&offering=community>
* Python 2.7 or above
    * If you don't have any version of python installed, Python 3.6 is preferable
    * <https://www.python.org/downloads/>
* cwltool
    * A command-line executor for CWL
    * This has to be installed using the command line
    * <https://github.com/common-workflow-language/cwltool#install>
* A text editor for code
    * If you don't already have a favourite, I recommend Atom
    * <https://atom.io/>

## Slides
[Workshop Slides](media/index.html) (use the arrow keys to navigate)

* [Part 1: Introduction](media/index.html#3)
* [Part 2: Tools](media/index.html#7)
* [Part 3: Writing Workflows](media/index.html#42)
* [Part 4: YAML](media/index.html#51)
