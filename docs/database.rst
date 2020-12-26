database
==================

.. py:module::src.database

This module acts on the Database namedtuple object. We do this C-esque work 
instead of a python class for the slight memory efficiency and speed bump. 

.. autofunction:: src.database.extract_protein_name

--------------------

.. autofunction:: src.database.build

--------------------

.. autofunction:: src.database.get_proteins_with_subsequence

--------------------

.. autofunction:: src.database.get_proteins_with_subsequence_ion

--------------------

.. autofunction:: src.database.get_entry_by_name