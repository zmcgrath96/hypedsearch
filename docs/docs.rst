docs
====

This page will walk you through what tools are used to make the docs, 
how to make the docs, and how to keep them updated

Tools
^^^^^
These docs are created using sphinx_ and ReadTheDocs_. While a little 
tricky to get the hang of, it allowed for this wiki page. If you want 
to learn how it was initially set up, follow along the video tutorial in the 
Read The Docs page. 

.. _sphinx: https://www.sphinx-doc.org/en/master/
.. _ReadTheDocs: https://readthedocs.org/

Building
^^^^^^^^
To build the docs, make sure you have sphinx_ installed. Follow their 
installation instructions, and if you have any problems, use this `stack overflow`__ 
question for troubleshooting. 

__ https://stackoverflow.com/questions/37757151/the-sphinx-build-command-was-not-found

To build the docs, do the following 

.. code-block:: bash

    $> cd docs
    $docs> make html

If there are not errors then you simply need to open the `docs/_build/html/index.html` 
file to see the docs locally. 

Updating
^^^^^^^^
These docs are hosted for free in ReadTheDocs, and a commit hook is used to 
update the build and update the logs to latest. Simply update the appropriate 
python docstrings and .rst files and push to see the updates.