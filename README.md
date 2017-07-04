

# Bioinformatics Tutorials and Workshops

These tutorials have been developed by Melbourne Bioinfomatics (formerly VLSCI) for training in commmon bioinformatics tasks.

The tutorials have been deployed here: http://vlsci.github.io/lscc_docs/tutorials/

Tutorials are written in [Markdown](http://en.wikipedia.org/wiki/Markdown) and built with `MkDocs`(http://www.mkdocs.org/).
All documentation source files are under `docs/`. Tutorials are stored under `docs/tutorials/`.

Slides for workshops are currently stored in Google Drive as Google Slides. Slides may optionally be added to this repository as e.g. Markdown or PDF.

## What's in a tutorial

Each tutorial subdirectory under `docs/tutorials/` should contain:

* At least one Markdown (`.md`) document with tutorial instructions. One of these should be the main tutorial document and be listed in `mkdocs.yml` in the top directory of this repository.
* A `README.md` Markdown document containing suggested PR reviewers and the location of any relevant Google Slides (or any other slides that cannot be stored in the repository directly).
* A `media` subdirectory with images that will be used by MkDocs.

Tutorials are built and deployed to the `gh-pages` branch using MkDocs, and then appear at http://vlsci.github.io/lscc_docs/tutorials/.

## How to contribute

### Fork and clone this repository

- Fork this repository in github.
- Clone your fork to your local computer, e.g.: `git clone https://github.com/<your_account>/lscc_docs`
- `cd lscc_docs`

### Sync with upstream

If the repository has been edited since you forked it, you will need to bring your fork up to date before you can contribute changes:

```
git remote add upstream https://github.com/vlsci/lscc_docs
git fetch upstream
git merge upstream/master
```

### Set up your environment and build the documentation

We are using the Python package [MkDocs](http://www.mkdocs.org/).

I recommend installing MkDocs into a virtualenv:

```
virtualenv mkdocs_dev
source mkdocs_dev/bin/activate
pip install -U mkdocs markdown-include
pip install -U mdx_showable
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

## Deploy to your fork (optional)

You can deploy changes to your fork. This will allow you (and PR reviewers) to
view your changes online prior to merging your pull request.

Run:

```
mkdocs gh-deploy
```

This will build the Markdown into HTML in your fork's `gh-pages` branch, AND immediately push the result without any further chance to review it. Your build should be
visible at `http://<youraccount>.github.io/lscc_docs/tutorials`.

TODO: check steps needed to have a github.io site up.

### Merging and deploying a pull request

Ideally you should have someone else check and merge your pull request rather than do it yourself. Suggested reviewers for each tutorial are in that tutorial's `README.md`.

To merge and deploy a PR:

* Check the changes in github. Optionally, view the deployed site at the requester's fork.
* Merge the PR in github. If there are merge conflicts, ask the issuer of the pull request to bring their fork up to date (see "Sync with upstream", above) and re-issue the pull request.
* Clone or update a local copy of this repository, and update deployment of documentation, with:

```
git clone https://github.com/vlsci/lscc_docs
cd lscc_docs
mkdocs gh-deploy
```

or if you already have a local clone:

```
git checkout master
git pull
mkdocs gh-deploy
```

Check that the updated tutorial appears under http://vlsci.github.io/lscc_docs/tutorials/.

## Making changes to tutorial instructions

To edit tutorial instructions, just edit the `.md` files containing those instructions and commit your changes with a helpful commit message. Follow the instructions above to contribute changes.

## Making changes to slides, or adding slide sets

For slides stored in Google Drive, a link should be recorded in `docs/tutorials/<tutorialname>/README.md`. This link should point to the Google slides source (not PDF) where possible.

If you have made an alternate version of the slides for a workshop, you can list the link for it under "Other slides". Give some kind of description for this alternate version, e.g. "slides for a clinical audience", "slides for bioinformaticians at GCC", or "Clare's abbreviated 5-minute slides".

If you want to update the latest set of slides for a workshop:

* Move the link from "Current slides" in `README.md` to a bullet point under "Other slides". Give some kind of description for this old version, even if it is just "Slides prior to June 2017 edits".
* Create the new set of slides. Ideally, it should be stored in the LSCC shared Google Drive folder, in an appropriately named subfolder under either  `LSCC_shared/capacity_building/tutorials_workshops/` or `LSCC_shared/capacity_building/LSCC_NGSschool/`. Add the link to your new slides after "Current slides:" in the `README.md`.

Don't be afraid to replace the current slides with your version if you think it is more up to date. The old slides are still accessible.

## Adding a new tutorial

To create your new tutorial:

* Create a new subdirectory under `docs/tutorials` with a meaningful name, e.g. `docs/tutorials/molecular_modelling`
* Create a main Markdown document within that folder, containing the tutorial instructions, e.g. `docs/tutorials/molecular_modelling/molecular_modelling.md`. You can use other tutorials as a template for the layout.
* Add your main Markdown document to `mkdocs.yml` in the correct section. This will cause a link to it to appear in the menu of tutorials.
* Create a `media` subdirectory, e.g. `docs/tutorials/molecular_modelling/media`. Copy in any images you need from other tutorials (e.g. logos). You can add any images here that you want to link into your Markdown documents.
* Create a `README.md` file, e.g. `docs/tutorials/molecular_modelling/README.md`. You can copy the format from another tutorial. Add your own name as a PR reviewer for future changes.

Add slides:

* If you are using Google Slides, ideally store them in a new subdirectory under `LSCC_shared/capacity_building/LSCC_NGSschool`.
* Add a link to the slides to your new `README.md` with the text "Current slides".

Once you've created your content:

* Follow the "How to contribute" instructions above to preview your changes locally. Check that the new tutorial appears in the menu and renders correctly. Then pull request and (ideally) ask someone to merge and deploy your changes.
* Tell Christina that there is a new workshop available.
