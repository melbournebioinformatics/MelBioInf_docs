# LSCC Documentation

## How to build the documentation

We are using the Python package [MkDocs](http://www.mkdocs.org/).

I recommend installing MkDocs into a virtualenv:

```
virtualenv mkdocs_dev
source mkdocs_dev/bin/activate
pip install -U mkdocs
pip install -U mdx_showable
```

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

You can build the site using the command:

```
mkdocs build
```

This will put all the HTML, CSS, Javascript etcetera for the site in the directory called *site*.

## How to edit the pages

The pages are written in [Markdown](http://en.wikipedia.org/wiki/Markdown).
