

# Bioinformatics Tutorials

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Bioinformatics Tutorials](#bioinformatics-tutorials)
	- [What's in a tutorial](#whats-in-a-tutorial)
	- [How to contribute](#how-to-contribute)
		- [Fork and clone this repository](#fork-and-clone-this-repository)
		- [Sync with upstream](#sync-with-upstream)
		- [Set up your environment and build the documentation](#set-up-your-environment-and-build-the-documentation)
		- [Deploy to your fork (optional)](#deploy-to-your-fork-optional)
		- [Pull request](#pull-request)
	- [Merging and deploying a pull request](#merging-and-deploying-a-pull-request)
	- [Making changes to tutorial instructions](#making-changes-to-tutorial-instructions)
	- [Making changes to slides](#making-changes-to-slides)
	- [Adding a new tutorial](#adding-a-new-tutorial)

<!-- /TOC -->

These tutorials have been developed by Melbourne Bioinformatics (formerly VLSCI) and are used in Melbourne Bioinformatics hands-on workshops.

The tutorials have been deployed here: http://melbournebioinformatics.github.io/MelBioInf_docs/

Tutorials are written in [Markdown](http://en.wikipedia.org/wiki/Markdown) and built with [MkDocs](http://www.mkdocs.org/).
All documentation source files are under `docs/`. Tutorials are stored under `docs/tutorials/`.

Slides for workshops are currently stored in Google Drive as Google Slides. Slides may optionally be added to this repository as e.g. Markdown or PDF.

## What's in a tutorial

Each tutorial subdirectory under `docs/tutorials/` should contain:

* At least one Markdown (`.md`) document with tutorial instructions. One of these should be the main tutorial document and be listed in `mkdocs.yml` in the top directory of this repository.
* A `README.md` Markdown document containing suggested PR reviewers and the location of any relevant Google Slides (or any other slides that cannot be stored in the repository directly).
* A `media` subdirectory with images that will be used by MkDocs.

Tutorials are built and deployed to the `gh-pages` branch using MkDocs, and then appear at http://melbournebioinformatics.github.io/MelBioInf_docs/.

## How to contribute

This section is a guide to contributing edits to the repository. For a guide to which files to edit, see [Making changes to tutorial instructions](#making-changes-to-tutorial-instructions), [Making changes to slides](#making-changes-to-slides), and [Adding a new tutorial](#adding-a-new-tutorial).

### Fork and clone this repository

You should do your work in a fork under your own github account:

* Fork this repository in github using the fork button at https://github.com/melbournebioinformatics/MelBioInf_docs.
* Clone your fork to your local computer, e.g.: `git clone https://github.com/<your_account>/MelBioInf_docs`

It is often a good idea to make a separate branch in your fork to work on your changes.

### Sync with upstream

If the repository has been edited since you forked it, you will need to bring your fork up to date before you can contribute changes:

```
git remote add upstream https://github.com/melbournebioinformatics/MelBioInf_docs
git fetch upstream
git merge upstream/master
```

### Set up your environment and build the documentation

We are using the Python package [MkDocs](http://www.mkdocs.org/).

I recommend installing MkDocs into a virtualenv:

```
python3 -m venv mkdocs_dev # or python2 -m virtualenv mkdocs_dev
source mkdocs_dev/bin/activate
pip install -r requirements.txt
```

Once you have mkdocs installed (and in your PATH) then you can preview the website by running the following command in the top directory of this repository:

```
mkdocs serve
```

This will start up a web server hosting on a local URL, like so:
```
Running at: http://127.0.0.1:8000/
Live reload enabled.
Hold ctrl+c to quit.
```

You can view the site if you point your browser at the specified URL.

MkDocs will automatically try to update the local site preview if you edit the documentation pages.

Alternatively, you can build the site locally using `mkdocs build`. This will put all the HTML, CSS, Javascript etcetera for the site in the directory called *site*, but unlike `mkdocs serve`, will not run a local webserver or update the build when you make changes.

It is a good idea to preview your changes locally before pushing them.

### Deploy to your fork (optional)

You can deploy changes to your fork. This will allow you (and PR reviewers) to
view your changes online prior to merging your pull request.

Run:

```
mkdocs gh-deploy
```

This will build the Markdown into HTML in your fork's `gh-pages` branch, AND immediately push the result without any further chance to review it. Your build should be
visible at `http://<your_account>.github.io/lscc_docs/tutorials`.


### Pull request

Commit your changes. Then push them to your fork with `git push origin master`.

Open your fork in github at `https://github.com/<your_account>/MelBioInf_docs` and create a pull request.

If your fork is up to date, you should see "These branches can be automatically merged" while you are creating the pull request. If your fork is *not* up to date you should bring it up to date (see [Sync with upstream](#sync-with-upstream))) and resolve any merge conflicts before creating the pull request.

Ideally you should have someone else check and merge your pull request rather than do it yourself. Suggested reviewers for each tutorial are in that tutorial's `README.md`.

Preferably, change one major thing per pull request - e.g. edit one tutorial in one pull request, and make a separate pull request if you want to edit another tutorial.

## Merging and deploying a pull request

If someone has asked you to merge their PR, or if you are merging in your own:

* Check the changes in github to spot any errors. You should be able to see the diff for the changes if the requester has sent you a link to their Pull Request.
* Optionally, to view the new docs in their final form:
	- view the deployed site at the requester's fork. This will be at `http://<requesters_account>.github.io/MelBioInf_docs/` *if* the author has run `mkdocs gh-deploy` on their fork. Or,
	- clone and build the requester's fork locally, using the [instructions above](#set-up-your-environment-and-build-the-documentation).
* Merge the PR in github. If there are merge conflicts, ask the issuer of the pull request to bring their fork up to date (see [Sync with upstream](#sync-with-upstream)) and re-issue the pull request.
* Clone or update a local copy of this repository, and re-deploy the updated documentation, like so:

```
git clone https://github.com/melbournebioinformatics/MelBioInf_docs
cd MelBioInf_docs
mkdocs gh-deploy
```

Check that the updated tutorial appears under http://melbournebioinformatics.github.io/MelBioInf_docs/

## Making changes to tutorial instructions

Tutorial instructions are stored as Markdown and fully versioned, so you can just edit the `.md` files containing those instructions and commit your changes. Follow the instructions above to contribute changes.

New media can be added to the tutorial's `media` subdirectory and linked in to the Markdown document.

## Making changes to slides

For slides stored in Google Drive, a link should be recorded in `docs/tutorials/<tutorialname>/README.md`. This link should point to the Google slides source (not PDF) where possible.

If you have made an alternate version of the slides for a workshop, you can list the link for it under "Other slides" in the `README.md`. Give some kind of description for this alternate version, e.g. "slides for a clinical audience", "slides for bioinformaticians at GCC", or "Clare's abbreviated 5-minute slides".

If you want to update the latest set of slides for a workshop:

* Move the link from "Current slides" in `README.md` to a bullet point under "Other slides". Give some kind of description for this old version, even if it is just something like "Slides prior to June 2017 edits".
* Create the new set of slides. Ideally, it should be stored in the LSCC shared Google Drive folder, in an appropriately named subfolder under either  `LSCC_shared/capacity_building/tutorials_workshops/` or `LSCC_shared/capacity_building/LSCC_NGSschool/`.
* Get a shareable link to your slides from Google Drive. If you do this via "Get shareable link", you will find that Google will automatically turn on link sharing. You should make sure to turn off edit permissions on link sharing or, if preferred, turn off link sharing completely: the link will still work, and so long as the slides were created in the `LSCC_shared` folder, they will still be accessible to anyone with access permissions for that folder.
* Add the link to your new slides after "Current slides:" in the `README.md`.

Don't be afraid to replace the current slides with your version if you think it is more up to date. The old slides are still accessible.

## Adding a new tutorial

To create your new tutorial:

* Create a new subdirectory under `docs/tutorials` with a meaningful name, e.g. `docs/tutorials/molecular_modelling`
* Create a main Markdown document within that folder, containing the tutorial instructions, e.g. `docs/tutorials/molecular_modelling/molecular_modelling.md`. You can use other tutorials as a template for the layout.
* Add your main Markdown document to `mkdocs.yml` in the correct section. This will cause a link to it to appear in the menu of tutorials.
* Create a `media` subdirectory, e.g. `docs/tutorials/molecular_modelling/media`. Copy in any images you need from other tutorials (e.g. logos). You can add any images here that you want to link into your Markdown documents.
* Create a `README.md` file in your tutorial subdirectory, e.g. `docs/tutorials/molecular_modelling/README.md`. You can copy the format of this file from another tutorial. Add your own name as a PR reviewer for future changes.

Add slides:

* If you are using Google Slides, ideally store them in a new subdirectory under `LSCC_shared/capacity_building/LSCC_NGSschool`.
* Add a link to the slides to your new `README.md` and label them "Current slides".

Once you've created your content:

* Before committing changes, add all newly created files to git with `git add`.
* Follow the [How to contribute](#how-to-contribute) instructions above to create a pull request. When previewing your changes, check that the new tutorial appears in the menu and renders correctly.
* Tell Christina that there is a new workshop available.
