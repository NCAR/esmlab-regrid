============================
Contributing to esmlab-regrid
============================

Contributions are highly welcomed and appreciated.  Every little help counts,
so do not hesitate!

.. contents:: Contribution links
   :depth: 2


.. _submitfeedback:

Feature requests and feedback
-----------------------------

Do you like esmlab-regrid?  Share some love on Twitter or in your blog posts!

We'd also like to hear about your propositions and suggestions.  Feel free to
`submit them as issues <https://github.com/NCAR/esmlab-regrid>`_ and:

* Explain in detail how they should work.
* Keep the scope as narrow as possible.  This will make it easier to implement.


.. _reportbugs:

Report bugs
-----------

Report bugs for esmlab-regrid in the `issue tracker <https://github.com/NCAR/esmlab-regrid>`_.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting,
  specifically the Python interpreter version, installed libraries, and esmlab-regrid
  version.
* Detailed steps to reproduce the bug.

If you can write a demonstration test that currently fails but should pass
(xfail), that is a very useful commit to make as well, even if you cannot
fix the bug itself.


.. _fixbugs:

Fix bugs
--------

Look through the `GitHub issues for bugs <https://github.com/NCAR/esmlab-regrid/labels/type:%20bug>`_.

Talk to developers to find out how you can fix specific bugs.


Write documentation
-------------------

esmlab-regrid could always use more documentation.  What exactly is needed?

* More complementary documentation.  Have you perhaps found something unclear?
* Docstrings.  There can never be too many of them.
* Blog posts, articles and such -- they're all very appreciated.

You can also edit documentation files directly in the GitHub web interface,
without using a local copy.  This can be convenient for small fixes.

.. note::
    Build the documentation locally with the following command:

    .. code:: bash

        $ conda env update -f ci/environment-dev-3.7.yml
        $ cd docs
        $ make html

    The built documentation should be available in the ``docs/_build/``.



 .. _`pull requests`:
.. _pull-requests:

Preparing Pull Requests
-----------------------


#. Fork the
   `esmlab-regrid GitHub repository <https://github.com/NCAR/esmlab-regrid>`__.  It's
   fine to use ``esmlab-regrid`` as your fork repository name because it will live
   under your user.

#. Clone your fork locally using `git <https://git-scm.com/>`_ and create a branch::

    $ git clone git@github.com:YOUR_GITHUB_USERNAME/esmlab-regrid.git
    $ cd esmlab-regrid

    # now, to fix a bug or add feature create your own branch off "master":

    $ git checkout -b your-bugfix-feature-branch-name master


#. Install `pre-commit <https://pre-commit.com>`_ and its hook on the esmlab-regrid repo::

     $ pip install --user pre-commit
     $ pre-commit install

   Afterwards ``pre-commit`` will run whenever you commit.

   https://pre-commit.com/ is a framework for managing and maintaining multi-language pre-commit hooks
   to ensure code-style and code formatting is consistent.
#. Install dependencies into a new conda environment::

    $ conda env update -f ci/environment-dev-3.7.yml


#. Run all the tests

   Now running tests is as simple as issuing this command::

    $ conda activate esmlab-regrid-dev
    $ pytest --junitxml=test-reports/junit.xml --cov=./ --verbose


   This command will run tests via the "pytest" tool against Python 3.7.

#. You can now edit your local working copy and run the tests again as necessary. Please follow PEP-8 for naming.

   When committing, ``pre-commit`` will re-format the files if necessary.

#. Commit and push once your tests pass and you are happy with your change(s)::

    $ git commit -a -m "<commit message>"
    $ git push -u


#. Finally, submit a pull request through the GitHub website using this data::

    head-fork: YOUR_GITHUB_USERNAME/esmlab-regrid
    compare: your-branch-name

    base-fork: NCAR/esmlab-regrid
    base: master
