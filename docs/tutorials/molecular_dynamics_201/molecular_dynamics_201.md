<p>
<img src="../../../img/melbioinf_logo.png" alt="Melbourne Bioinformatics logo" align="left" width="40%"/><br />
</p>
<p><br></p>
<p><br></p>

# Molecular Dynamics Tutorial

-----

Welcome to the practical session of the supercomputer cluster workshop. In the following exercises we will be logging on a high performance computer (HPC) submitting NAMD\* molecular dynamics jobs. As students have a vast range of skill levels and requirements, each of the following exercises are designed to be stand alone and worked on at your own pace. Start at where you feel most comfortable and skip the exercises that are too simple. As we will be working from a terminal on the cluster and later downloading data back to our local desktop for visualization and analysis, we will be assuming that the users have basic knowledge of Linux and the molecular visualization program VMD\*.

As we will be using a computer terminal for these exercises, it will be useful if you are familiar with some Unix/Linux commands. If you are unfamiliar with these, we suggest you first work through the exercises outlined on this [Unix tutorial page](https://melbournebioinformatics.github.io/MelBioInf_docs/tutorials/unix/robinson-unix-link/).

\**NAMD - molecular dynamics program*

\**VMD - molecular visualization program*

-----

## Exercise overview

1. [Basic introduction to cluster computing]()

If you have never launched a job on a cluster before, this exercise is for you. Here we will simply copy across a bare minimum directory containing input files ready to go for a short NAMD simulation. Once the job is finished, we will download the data to your local computer and visualize the trajectory with VMD. (this exercise starts in this document)

2. Intermediate molecular dynamics exercise - building input files

In this exercise we will be introducing a more sophisticated directory structure using scripts to direct the output and run the job. We initially download a copy of the directory to our local computer where we then build up the input files using VMD. We then upload this directory to the cluster where we submit the job. Finally we will download the results back to our local computer to visualize with VMD. (This exercise starts in the document: Intermediate_MD_workshop)
