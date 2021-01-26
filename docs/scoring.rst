scoring
==============

.. py:module::src.scoring

The two modules initially were suppossed to be different, but since have become
nearly one in the same. Compare masses is the real worker for scoring spectra, 
and scoring is an interface to that as well as other scoring metrics. Mass 
comparisons is still separate due to all of the experimentation tried to 
make the scoring as fast as possible. For both files, all older attempts have 
just been commented out for reference so that we don't try the same things over 
again. 


mass_comparisons
^^^^^^^^^^^^^^^^

.. autofunction:: src.scoring.mass_comparisons.optimized_compare_masses


scoring
^^^^^^^

.. autofunction:: src.scoring.scoring.score_sequence

--------------------

.. autofunction:: src.scoring.scoring.hybrid_score

--------------------

.. autofunction:: src.scoring.scoring.precursor_distance

--------------------

.. autofunction:: src.scoring.scoring.total_mass_error

--------------------

.. autofunction:: src.scoring.scoring.digest_score
