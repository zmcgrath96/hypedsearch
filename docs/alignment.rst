alignment
=========

.. py:module::src.alignment

This is the heaviest lifter of the entire project. The alignment process happens
here. The input is a spectrum, and the output is a set of alignments.

alignment_utils
^^^^^^^^^^^^^^^

.. autofunction:: src.alignment.alignment_utils.__split_hybrid

--------------------

.. autofunction:: src.alignment.alignment_utils.__get_surrounding_amino_acids

--------------------

.. autofunction:: src.alignment.alignment_utils.__add_amino_acids

--------------------

.. autofunction:: src.alignment.alignment_utils.__remove_amino_acids

--------------------

.. autofunction:: src.alignment.alignment_utils.align_overlaps

--------------------

.. autofunction:: src.alignment.alignment_utils.match_precursor

--------------------

.. autofunction:: src.alignment.alignment_utils.get_parents

--------------------

.. autofunction:: src.alignment.alignment_utils.extend_non_hybrid 


alignment
^^^^^^^^^

.. autofunction:: src.alignment.alignment.same_protein_alignment

--------------------

.. autofunction:: src.alignment.alignment.align_b_y

--------------------

.. autofunction:: src.alignment.alignment.attempt_alignment


hybrid_alignment
^^^^^^^^^^^^^^^^

.. autofunction:: src.alignment.hybrid_alignment.__replace_ambiguous_hybrid

--------------------

.. autofunction:: src.alignment.hybrid_alignment.replace_ambiguous_hybrids

--------------------

.. autofunction:: src.alignment.hybrid_alignment.hybrid_alignment