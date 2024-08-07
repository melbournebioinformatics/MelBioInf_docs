# Version control with Git

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
*Section 1:* What is Git? And why use it?  
*Section 2:* Creating your first Git repository  
*Section 3:* Incorporating Git into your day-to-day work  
*Section 4:* Leveraging the best out of GitHub

**Learning objectives:** learn the basics of version control, how to track changes and how to collaborate.

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

??? note "Configuring SSH access for GitHub"
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

### 1. What is Git?
A concise definition is that Git is a **version control system**. It was developed in the early 2000s by 
Linus Torvalds, the person behind the Linux kernel, and other members of the Linux community. In fact, the initial purpose of
Git was to support the development and updates of the Linux kernel.

Git is attractive because of several things:

* It's **fully distributed**. What this means is that several (hundreds, thousands) of users can work on the same project in a streamlined way.
* It's **blazing fast**. Because (nearly) every operation that Git performs is local, everything works almost instantaneously.
* It has **integrity**. Git checksums data before storing it, so it's impossible to change the content of files without knowing about it.

This is all well and good, but let's try to understand what version control means for us as researchers.

![](http://www.phdcomics.com/comics/archive/phd101212s.gif)

### 2. Why use version control?
- Backing up your code
- Sharing your code
    - Packaging and distribution
- Documenting your work
- Collaborating
    - The "Lingua Franca" of software engineering
    - Will be required if you want to work in coding projects with other people (including past and future you!)
- Why Git will make you a better researcher
    - You know you can always go back to a working version of your code
    - You will have a digital lab notebook 
    - You will have a way of showcasing your projects
    - You will be able to distribute your code to others
    - You will be able to modify other peoples' code and make contributions to it

### 3. Basic concepts
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