# Bioinformatics Documentation

These tutorials have been developed by Melbourne Bioinfomatics (formerly VLSCI) for training in commmon bioinformatic tasks.

## Deployment

The tutorials have been deployed here: http://vlsci.github.io/lscc_docs/tutorials/

## How to contribute

## Install the mkdocs tools

We are using the Python package [MkDocs](http://www.mkdocs.org/).

I recommend installing MkDocs into a virtualenv:

```
virtualenv mkdocs_dev
source mkdocs_dev/bin/activate
pip install -U mkdocs markdown-include
pip install -U mdx_showable
```

## Fork and clone the repo

- Fork this repository to your own GitHub repository
- Git clone to your local computer: `git clone <your repo name> melbio_tutorials`
- `cd melbio_tutorials`
- Make changes, e.g. in Atom, a Markdown editor

The pages are written in [Markdown](http://en.wikipedia.org/wiki/Markdown).

- Git add, commit, push

Once you have mkdocs installed (and in your PATH) then you can preview the website by running the
following command in the top directory of this repository:

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

MkDocs will automatically try to update the site if you edit the documentation pages.

- Within GitHub, send a pull request to MB
- Suggested Pull Request reviewers to contact for each section

## Sync with upstream

- Point to the upstream repository: `git remote add upstream https://github.com/vlsci/lscc_docs.git `
- Fetch upstream changes: `git fetch upstream`
- Merge with your copy: `git merge upstream/master`

## Add a new topic

The pages are in the docs/tutorials folder.

To add a new topic, create a folder within docs/tutorials. e.g. "Assembly"

Create a page within that folder, e.g. Assembly/index.md

Add the file to the mkdcos.yml in the correct section.

## Deploy to the public website


You can build the site using the command:

```
mkdocs build
```

This will put all the HTML, CSS, Javascript etcetera for the site in the directory called *site*.


## Slides

- Where are they located - google drive folder

    - LSCC_shared/capacity_building/

- Main and alternate sets
- How to update the slides / add alternate slides
