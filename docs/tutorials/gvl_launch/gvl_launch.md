<style type="text/css">
    body{
        line-height: 2;
        font-size: 16px;
    }

    ol li{padding: 4px;}
    ul li{padding: 0px;}
    h4 {margin: 30px 0px 15px 0px;}

    div.code {
        font-family: "Courier New";
        border: 1px solid;
        border-color: #999999;
        background-color: #eeeeee;
        padding: 5px 10px;
        margin: 10px;
        border-radius: 5px;
        overflow: auto;
    }

    div.question {
        color: #666666;
        background-color: #e1eaf9;
        padding: 15px 25px;
        margin: 20px;
        font-size: 15px;
        border-radius: 20px;
    }

    div.question h4 {
        font-style: italic;
        margin: 10px 0px !important;
    }
</style>

<img src="../media/gvl_logo.jpg" height=100px align=right>

# Launching a Personal GVL Server on the NeCTAR Research Cloud

-----

## Tutorial Overview

This document guides you to launch your own GVL Analysis platform, with Galaxy,
using your default NeCTAR allocation.

The tutorial will go over:

1. Accessing the NeCTAR dashboard using Australian Access Federation (AAF)
   credentials
2. Getting your NeCTAR Research Cloud credentials
3. Launching the GVL image
4. Accessing your GVL instance
5. GVL services
6. Shutting your machine down

-----

## Background

#### What is the GVL?
The Genomics Virtual Laboratory (GVL) is a research platform that can be
deployed on the cloud.

<img src="../media/gvl_dashboard.png" height=500px style="display:block; margin-left: auto; margin-right:auto;">

A private GVL server is a virtual machine running on the cloud and contains a
pre-installed suite of tools for performing bioinformatics analyses. It differs
from public GVL servers (such as the
[Galaxy Tutorial Server](http://galaxy-tut.genome.edu.au/),
[Galaxy Melbourne](http://galaxy-mel.genome.edu.au/), and
[Galaxy Queensland](http://galaxy-qld.genome.edu.au/)) by providing full
administrative access to the server, as well as the full suite of GVL services,
whereas public GVL servers provide restricted access for security reasons.
For example, public GVL servers do not provide access to the Ubuntu desktop,
the Linux command line or JupyterHub at present.

Accessing the GVL server is completely free on the Australian NeCTAR Research
Cloud, provided that you have a account with NeCTAR with allocated resources.

#### What is NeCTAR?

The [National eResearch Collaboration Tools and Resources project](http://nectar.org.au/)
(NeCTAR) is an Australian programme that provides computing infrastructure
and services to Australian researchers. The NeCTAR Cloud allows us to
deploy virtual machines as a platform for research.

While it is possible to launch the GVL on [Amazon](https://aws.amazon.com/),
you may have to pay Amazon usage charges (the GVL software itself is free).

-----

## Section 1: Access the NeCTAR dashboard

1.  Login into the NeCTAR dashboard at
    [dashboard.rc.nectar.org.au](https://dashboard.rc.nectar.org.au).  
    Only members of the Australian Access Federation (AAF) can access Australian
    Research Cloud resources. Most Australian universities are members of the
    AAF, however, if you belong to an institution that is not a member of AAF,
    [GVL Help](http://genome.edu.au/help/gvl-help) may be able to provide you
    with credentials. You also have the option of using a commercial cloud
    provider such as [AWS](https://aws.amazon.com/) to host your server.

2.  Choose your organisation from the list and login using your credentials.
    If you are from the University of Melbourne, choose 'The University of
    Melbourne' and not 'The University of Melbourne (with ECP)'.

3.  Login with your institutional username and password. If this is the first
    time you have accessed the Australian Research Cloud, you must agree to some
    terms and conditions.

    When you log in for the first time, you are automatically allocated a trial
    project which lasts for 3 months. This trial project allows you to launch a
    medium instance (with 2 cores) which is sufficient for launching a GVL
    instance.

    <img src="../media/nectar_dashboard.png" height=400px style="display:block; margin-left: auto; margin-right:auto;">

    If you have projects which require more compute resources, you can apply
    for more allocation [here](https://dashboard.rc.nectar.org.au/allocation/).

-----

## Section 2: Get your cloud credentials

Launching a GVL instance requires EC2 API keys from NeCTAR. Obtaining these keys
is a simple 3 step process.

<img src="../media/get_ec2_annotated.png" height=400px/>

1.  From the NeCTAR dashboard, on the left sidebar, navigate to  
    Project > Compute > Access & Security

2.  Click on the top **'API Access'** tab.

3.  Click on the **'+ View Credentials'** button on the top right. A window
    containing your API access key (1) and your secret key (2) will appear.
    Click on the eye icon to view your secret key. Keep this key secret and
    secure.

    <img src="../media/ec2_creds_annotated.png" height=300px/>

    Keep this page open for later reference.

-----

## Section 3: Launch your personal GVL instance

1.  In a new browser tab, go to [launch.genome.edu.au](https://launch.genome.edu.au)

2.  Fill in the required fields:
    - **Cloud:** Keep the default: NeCTAR (OpenStack)
    - **Access key / Secret key:** Copy and paste your access key and your secret
      key you obtained in the previous step in the required fields.
    - **Institutional email:** Enter your institutional email address
    - **Cluster name:** Choose 'Specify a new name' and enter a name for your
      instance (eg. GVL_workshop). It is recommended you choose a unique name
      if you launch multiple instances.
    - **Password:** Choose a strong password and remember it. This is the
      password you will use later to log into your instance.
    - **Instance type:** Keep the default: Medium (2 vcpu / 8GB RAM)
    - **Cluster type:** Keep the default: Cluster with Galaxy
    - **Storage type:** Keep the default: Transient instance storage

    <img src="../media/launch_gvl.png" height=500px style="display:block; margin-left: auto; margin-right:auto; margin-top:10px;">

3.  Optional advanced options  
    Toggle the 'Show advanced startup options' option to see more options. For
    this tutorial, it is not necessary to modify any of the advanced options.
    - **Placement:** You can also choose the region of where your server is
      hosted. If you are doing lots of data transfer, it may be beneficial to
      pick a location close to your physical location. More information about
      zones can be found [here](https://support.rc.nectar.org.au/docs/availability-zones).
    - **Key pair:** Key pairs are SSH credentials that can be used to access
      your instance. You can create and import key pairs in the NeCTAR dashboard
      by navigating to Project > Compute > Access & Security > Key Pairs and
      creating or importing a key pair.
    - **Image:** The image to start. The default option is usually the latest,
      most stable GVL image.
    - **Flavor:** Flavours are versions of the GVL with slightly customised
      toolsets. These toolsets are optimised for different usage patterns. More
      information about the different flavours can be found
      [here](http://genome.edu.au/get/get#gvlflavours). For this tutorial, keep
      the default option of 'GVL + Tutorial Indices'

    <img src="../media/launch_gvl_advanced_options.png" height=500px style="display:block; margin-left: auto; margin-right:auto; margin-top:10px;">

4.  Click **'Create a cluster'** to launch the image.  
    The launch process takes 2-5 minutes to start the machine and another 5
    minutes to start and configure Galaxy.

    <img src="../media/launch_gvl_wait.png" height=250px style="display:block; margin-left: auto; margin-right:auto;">

    If your progress bar seems stuck on the 'Requesting' stage for > 5 minutes,
    navigate back to the launch page and try selecting a different availability
    zone under the advanced startup options **Placement** field. You will need
    re-enter your secret key and your password.

-----

## Section 4: Access your GVL instance

1.  Once your instance has finished launching, click on the cluster IP address
    to access your GVL dashboard.

    <img src="../media/launch_gvl_ready.png" height=500px style="display:block; margin-left: auto; margin-right:auto;">

    If you accidentally closed the launch page, you can find your cluster's
    address on the [NeCTAR dashboard](https://dashboard.rc.nectar.org.au/project/instances/)
    by navigating to Project > Compute > Instances on the left panel. This page
    contains a list of your instances and can be used to terminate your instance
    if anything goes wrong. Copy and paste the instance's IP address into your
    browser's URL navigation bar.

2.  Explore the GVL dashboard.
    Have a read through of the services provided by the GVL.

<img src="../media/gvl_dashboard.png" height=500px style="display:block; margin-left: auto; margin-right:auto;">


-----

## Section 5: GVL services

Listed below are short descriptions of the services the GVL provides.

#### Galaxy
[Galaxy](https://galaxyproject.org/) is a web-based platform for computational
biomedical research. The GVL has a number of Galaxy tutorials available
[here](http://genome.edu.au/learn).

<img src="../media/gvl_galaxy.png" height=400px style="display:block; margin-left: auto; margin-right:auto;">

To begin using Galaxy, register as a user by navigating to User > Register
on the top Galaxy bar.

#### CloudMan
CloudMan is a cloud manager that manages services on your GVL instance. Use
Cloudman to start and manage your Galaxy service and to add additional nodes
to your compute cluster (if you have enough resources).

You can log into CloudMan by using the username 'ubuntu' and your cluster
password.

<img src="../media/cloudman.png" height=350px style="display:block; margin-left: auto; margin-right:auto;">

You can also shut down your instance (permanently) with CloudMan.

#### Lubuntu Desktop

Lubuntu is a lightweight desktop environment through which you can run desktop
applications on your virtual machine through your web browser. You can also
access the GVL command line utilities through the desktop.

You can log into Lubuntu Desktop using the username 'ubuntu' and your cluster
password.

<img src="../media/lubuntu.png" height=350px style="display:block; margin-left: auto; margin-right:auto;">

#### SSH

Secure Shell (SSH) is a network protocol that allows us to connect to a remotely
machine. You can login to your virtual machine remotely through an SSH client.

If you are using Windows, you will need to download an SSH client such as
[PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/). If you are using
OSX, open up a Terminal window.

If you are unfamiliar with the command line and UNIX, many
[tutorials](http://swcarpentry.github.io/shell-novice/)
[on UNIX](../unix/)
[can be](http://www.ee.surrey.ac.uk/Teaching/Unix/)
[found online](http://www.doc.ic.ac.uk/~wjk/UnixIntro/).

You can ssh into your machine using the either the username 'ubuntu' or the
username 'researcher' and using your cluster password. It is recommended to use
the researcher account when you are doing your computational research and use
the ubuntu account when you need administrative powers (such as installing
software).

<img src="../media/command_line.png" height=350px style="display:block; margin-left: auto; margin-right:auto;">

#### JupyterHub
[JupyterHub](http://jupyter.org/) is a web-based interactive computational
environment where you can combine code execution, text, mathematics, plots and
rich media into a single document. Currently, JupyterHub can connect to Python2
and Python3 kernals.

If you are unfamiliar with Python, there are many
[tutorials](http://swcarpentry.github.io/python-novice-inflammation/)
[available](../python_overview/python_overview/)
[online](http://pythonprogramminglanguage.com/).

You can log into JupyterHub with the username 'researcher' and your cluster
password.

You may need to install Python packages you intend to use via the command line
beforehand.

<img src="../media/jupyterhub.png" height=350px style="display:block; margin-left: auto; margin-right:auto;">

#### RStudio
RStudio Server gives browser-based access to RStudio, the popular programming
and analysis environment for the R programming language. You can find out more
about RStudio [here](https://www.rstudio.com/home/), and the R programming
language [here](https://www.r-project.org/).

You can log into RStudio with the username 'researcher' and your cluster
password.

<img src="../media/rstudio.png" height=350px style="display:block; margin-left: auto; margin-right:auto;">

#### Public HTML

This is a shared web-accessible folder. Any files you place in the directory
`/home/researcher/public_html` will be publicly accessible.

#### PacBio SMRT Portal
[PacBio's SMRT Portal](http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/) is an open source software suite for the analysis of single molecule, real-time sequencing data.

Before you use SMRT Portal, you need to firstly install it through the Admin console. *Please note PacBio recommends the use of a 16 core instance with 64GB of RAM (or higher) for this package.*

To install SMRT Portal from the GVL dashboard:

* Click on 'Admin' in the top navigation bar.

* Log in with the username 'ubuntu' and your cluster password. The screen will now show all the tools available to be installed.

* Scroll down to 'SMRT Analysis', and click 'install'.

* The install is complete when SMRT Portal is available as a tool on the GVL dashboard, and a green tick is displayed.

When launching SMRT Portal for the first time, you will need to register yourself as a new user.

-----

## Section 6: Shutting your machine down

There are two ways to terminate your instance. Terminating your instance is
permanent and all data will be deleted (unless you have persistent
[volume storage](https://support.rc.nectar.org.au/docs/volumes) which you will
need to apply for).

1.  Via Cloudman  
    Log into CloudMan by using the username 'ubuntu' and your cluster password.
    Click the **Shut down...** button under Cluster Controls.

2.  Via the NeCTAR dashboard  
    Navigate to the [Instances](https://dashboard.rc.nectar.org.au/project/instances/)
    page by navigating to Project > Compute > Instances on the left panel. Find
    the instance you want to terminate, and on the right-most column (Actions),
    click on the arrow button, and select **Terminate Instance**.
