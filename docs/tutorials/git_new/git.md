# Version Control with Git

## Overview

**Topic**

* [ ] Genomics
* [ ] Transcriptomics
* [ ] Proteomics
* [ ] Metabolomics
* [ ] Statistics and visualisation
* [ ] Structural Modelling
* [x] Basic skills


**Skill level**

* [x] Beginner  
* [ ] Intermediate  
* [ ] Advanced  

**Data:** None required, we will create our own code repository.

**Tools:** Git `>=2.40`

**Pipeline:**  
*Section 1:* Introduction  
*Section 2:* Your first Git repository  
*Section 3:* Collaborating with Git  
*Section 4:* Tips and best practices  

**Learning objectives:** learn the basics of version control, how to track changes in a repository and how to collaborate using Git.

!!! warning "Disclaimer"
    This tutorial is partially based on the [Software Carpentry's Git Novice lesson](https://swcarpentry.github.io/git-novice/). Some portions have been directly copied from that lesson, while others have been fully rewritten.
    The lesson is under a [CC-BY 4.0 license](https://swcarpentry.github.io/git-novice/LICENSE.html).

## Setup
**Install Git**  

If you have Anaconda, changes are you already have Git installed. Open either a Terminal or an Anaconda 
Prompt, and type `git --version`. If you are already have Git installed, it will show the version number.
Otherwise, you can install it with `conda install git -y`. If you don't have Anaconda or conda installed,
you can follow the instructions for [your operating system here.](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

**Create a GitHub account**  

The next step will be to create a GitHub account. If you already have one, skip this step. Otherwise, head
over to [GitHub](https://github.com/) and click the "Sign Up" button at the top right. Make sure you [verify
your email address.](https://docs.github.com/en/account-and-profile/setting-up-and-managing-your-personal-account-on-github/managing-email-preferences/verifying-your-email-address#troubleshooting-email-verification)

!!! success "You are all set for now."
    The two steps above are all that are required **before** the workshop. You can complete the following steps during the workshop.

**SSH access (optional)**

We strongly recommend that *at some point* you [configure SSH access to your GitHub account.](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
Otherwise, it will ask for your password every time you try to interact with GitHub. This is OK, but can 
be annoying after some time. See below on how to connect with SSH.

??? info "Configuring SSH access for GitHub"
    1. Make sure you have `openssl` installed by typing `openssl --version` in the Terminal.
        - You can install it with `conda install openssl` if not installed.
    2. Make sure your ssh-agent is running with ```eval `ssh-agent` ```.
    3. Create a new SSH key with `ssh-keygen`.
    4. Add the key to your SSH-agent with `ssh-add <PATH-TO-KEY>`, for example `ssh-add /home/username/.ssh/id_rsa`.
    5. Print your **public** key by typing `cat <PATH-TO-KEY>.pub`, for example `cat /home/username/.ssh/id_rsa.pub`.
    6. On GitHub, click your profile icon on the top right, then Settings > SSH and GPG keys > New SSH key > paste your **public** key and give it a meaningful name.
    7. Test SSH access with `ssh -T git@github.com`. It should print a message like this: `Hi <USERNAME>! You've successfully authenticated, but GitHub does not provide shell access.`

**Configuring Git**  

The final step in setting up Git is adding your user name and email to the Git configuration. This is required
because Git documents our work, that is, who did what. It is also a requirement for your GitHub authentication,
especially if you haven't set up SSH access.

!!! warning "Important"
    Use the same email that you used to sign up to GitHub.

    If you prefer to keep your email private, you can [set up a noreply address with GitHub](https://docs.github.com/en/account-and-profile/setting-up-and-managing-your-personal-account-on-github/managing-email-preferences/setting-your-commit-email-address).

Go to your Terminal and type:

``` bash
git config --global user.name "FirstName LastName"
git config --global user.email "your@email.com"
```

!!! success "Well done!"
    You are ready to start using Git.

## Introduction

### What is Git?
A concise definition is that Git is a **version control system**. It was developed in the early 2000s by 
Linus Torvalds, the person behind the Linux kernel, and other members of the Linux community. In fact, the initial purpose of
Git was to support the development and updates of the Linux kernel.

Git is attractive because of several things:

* It is **fully distributed**. What this means is that several (hundreds, thousands) of users can work on the same project in a streamlined way.
* It is **blazing fast**. Because (nearly) every operation that Git performs is local, everything works almost instantaneously.
* It has **integrity**. Git checksums data before storing it, so it's impossible to change the content of files without knowing about it.

This is all well and good, but let's try to understand what version control means for us as researchers.

![](media/phd101212s.png)

“notFinal.doc” by Jorge Cham, [https://www.phdcomics.com](http://www.phdcomics.com/comics/archive/phd101212s.gif)

---

We've all been in this situation before: it seems unnecessary to have
multiple nearly-identical versions of the same document. Some word
processors let us deal with this a little better, such as Microsoft
Word's
[Track Changes](https://support.office.com/en-us/article/Track-changes-in-Word-197ba630-0f5f-4a8e-9a77-3712475e806a),
Google Docs' [version history](https://support.google.com/docs/answer/190843?hl=en), or
LibreOffice's [Recording and Displaying Changes](https://help.libreoffice.org/Common/Recording_and_Displaying_Changes).

Version control systems start with a base version of the document and
then record changes you make each step of the way. You can
think of it as a recording of your progress: you can rewind to start at the base
document and play back each change you made, eventually arriving at your
more recent version.

![](media/play-changes.svg){alt='A diagram demonstrating how a single document grows as the result of sequential changes'}

Once you think of changes as separate from the document itself, you
can then think about "playing back" different sets of changes on the base document, ultimately
resulting in different versions of that document. For example, two users can make independent
sets of changes on the same document.

![](media/versions.svg){alt='A diagram with one source document that has been modified in two different ways to produce two different versions of the document'}

Unless multiple users make changes to the same section of the document - a 
[conflict](./git.md#glossary) - you can
incorporate two sets of changes into the same base document.

![](media/merge.svg){alt='A diagram that shows the merging of two different document versions into one document that contains all of the changes from both versions'}

A version control system is a tool that keeps track of these changes for us,
effectively creating different versions of our files. It allows us to decide
which changes will be made to the next version (each record of these changes is
called a [commit](./git.md#glossary), and keeps useful metadata
about them. The complete history of commits for a particular project and their
metadata make up a [repository](./git.md#glossary).
Repositories can be kept in sync across different computers, facilitating
collaboration among different people.

??? quote "The Long History of Version Control Systems"

    Automated version control systems are nothing new.
    Tools like [RCS](https://en.wikipedia.org/wiki/Revision_Control_System), [CVS](https://en.wikipedia.org/wiki/Concurrent_Versions_System), or [Subversion](https://en.wikipedia.org/wiki/Apache_Subversion) have been around since the early 1980s and are used by
    many large companies.
    However, many of these are now considered legacy systems (i.e., outdated) due to various
    limitations in their capabilities.
    More modern systems, such as Git and [Mercurial](https://swcarpentry.github.io/hg-novice/),
    are *distributed*, meaning that they do not need a centralized server to host the repository.
    These modern systems also include powerful merging tools that make it possible for
    multiple authors to work on
    the same files concurrently.

!!! question "Paper Writing"

    - Imagine you drafted an excellent paragraph for a paper you are writing, but later ruin it. How would you retrieve the *excellent* version of your conclusion? Is it even possible?

    - Imagine you have 5 co-authors. How would you manage the changes and comments they make to your paper? If you use LibreOffice Writer or Microsoft Word, what happens if you accept changes made using the `Track Changes` option? Do you have a history of those changes?

    ??? example "Solution"

        - Recovering the excellent version is only possible if you created a copy
        of the old version of the paper. The danger of losing good versions
        often leads to the problematic workflow illustrated in the PhD Comics
        cartoon at the top of this page.
        
        - Collaborative writing with traditional word processors is cumbersome.
        Either every collaborator has to work on a document sequentially
        (slowing down the process of writing), or you have to send out a
        version to all collaborators and manually merge their comments into
        your document. The 'track changes' or 'record changes' option can
        highlight changes for you and simplifies merging, but as soon as you
        accept changes you will lose their history. You will then no longer
        know who suggested that change, why it was suggested, or when it was
        merged into the rest of the document. Even online word processors like
        Google Docs or Microsoft Office Online do not fully resolve these
        problems.

!!! note "Keypoints"

    - Version control is like an unlimited ‘undo’.
    - Version control also allows many people to work in parallel.

### Why use version control?

Now that we've learned what version control **is**, let's understand why we should use it as researchers.

There are a number of reasons to argue that using version control will make us better researchers. That is because 
version control systems can be used for:

- **Backing up your code**  
The most fundamental idea of version control is that you can use it to safely back up your code. That means not only
having a copy of it, but having a copy of each version as your code evolves throughout time. If you use version control
effectively, it will be very difficult for your code to be permanently lost, deleted or erased.

- **Sharing your code**  
Although Git is already explicitly designed to work in a distributed manner, modern version control platforms
make sharing code with others even easier. If our repository is public, anyone can easily access it, copy, and modify the code as they please.
You can make your GitHub profile a portfolio of your coding and analysis projects, with a user-friendly interface.

    !!! example "Packaging and distribution"

        Platforms such as GitHub provide a number of features that facilitate **packaging** our code. One thing is to have a bunch of scripts in
        a repository, but if we want to distribute our code effectively, to make it easier for users to acquire and install our code, we can bundle
        it as a **software package** and upload it to platforms such as [PyPI](https://pypi.python.org/) (for Python) or [CRAN](https://cran.r-project.org/) (for R).

- **Collaborating**  
Sharing code with others is one thing, but Git also enables researchers to work *together* on the same project. You can review other users' commits,
and selectively apply or reject changes that they propose. This is further enabled by GitHub, which makes it easy to do it through the web interface.
Moreover, you can create organisations to host multiple projects, give collaborators write and admin access to projects, and give access to private repositories.

    !!! tip "The 'Lingua Franca' of software engineering"
        
        Because of the way Git enables collaboration, it has essentially become the way that programmers interact on a technical basis.
        If you want to make a contribution to a large code base or project, most likely you will have to submit your changes through Git.
        GitHub also allows the creation of **issues**, where you can report problems or create discussions about the code. Knowing Git will
        probably be required if you want to work in coding projects with other people (including past and future you!)

- **Documenting your work**  
Because changes in Git are structured through commits, it is very straightforward to document our work as we go. Whenever we create a commit,
we must write a message that's attached to it (we'll learn more about that), which almost mandates that we document what we are doing. This will
create a [history](./git.md#glossary) of our work which can effectively be used as a **digital research notebook** if done correctly. The way Git works also allows the creation of [branches](./git.md#glossary) and [tags](./git.md#glossary), which can be used to keep track of different parts of the development. This is especially useful for large projects where many people may be working on different things in parallel.

!!! note "Keypoints"
    **Why Git will make you a better researcher:**
    
    - You know you can always go back to a working version of your code
    - You will have a way of showcasing your projects
    - You will be able to distribute your code to others
    - You will be able to modify other peoples' code and make contributions to it
    - You will have a digital lab notebook 

### Basic concepts
- Repositories
    - Remotes (basic-level only)
- Commits
    - Snapshots of a repository
- The staging area

## Your first Git repository

### The Git directory

### Tracking changes
- The staging area
    - `git add` and `git commit`
- What makes a good commit?
    - Writing good commit messages

### Remotes: pushing and pulling

## Collaborating with Git

### Cloning a repository

### Cloning versus forking

### Submitting a pull request

### Reviewing a pull request

## Tips and best practices

### How to incorporate Git in your day-to-day work

### Leveraging the GitHub interface

### Do's and Dont's

#### Do
- Write good commit messages
    - Say **why**, not **what** you changed
- Make small changes
- Commit often
- Think about the reviewer
- Document as you code

#### Don't
- Make vague commit messages
- Accumulate unrelated changes in a single commit
- Let things go stale – delete or "stash" them

## Glossary

**Conflict:** A situation where changes from different commits or branches cannot be merged automatically.

**Commit:** A snapshot of changes to a repository, saved with a message describing the updates.

**Repository:** A storage space where your project’s files and version history are kept.

**History:** The record of all commits made in a repository, showing how the codebase evolved.

**Branch:** An independent line of development within a repository, allowing you to work on different features or fixes without affecting the main project.

**Tags:** Named markers pointing to specific commits, often used to label release versions.
