preprocessing
=============

.. py:module::src.preprocessing

The preprocessing module covers spectra loading and filtering, reducing the
search space over the database by tagging obseved masses with the k-mers with 
either a b or y ion mass that matches. The merge search is a linear time search 
for doing this matching process. As of now, this is not incredibly well optimized, 
and the limitation seems to be approximately 300 proteins to fit in 32 GB RAM 
comfortably. Future endeavors on preprocessing should try to reduce the time and 
space needed to tag these k-mers, possibly forgoing this approach altogether.

spectra_filtering 
^^^^^^^^^^^^^^^^^

.. autofunction:: src.preprocessing.spectra_filtering.relative_abundance_filtering

--------------------

.. autofunction:: src.preprocessing.spectra_filtering.peak_filtering

merge_search
^^^^^^^^^^^^

.. autofunction:: src.preprocessing.merge_search.merge

--------------------

.. autofunction:: src.preprocessing.merge_search.make_database_set

--------------------

.. autofunction:: src.preprocessing.merge_search.match_masses

preprocessing_utils
^^^^^^^^^^^^^^^^^^^

.. autofunction:: src.preprocessing.preprocessing_utils.load_spectra
