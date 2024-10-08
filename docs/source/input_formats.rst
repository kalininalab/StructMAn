Input Formatting
================

.. _smlf:

Simple Mutation List Format
----------------

Bla

.. _fasta:

Fasta Format
------------

| StructMAn can process a slightly modified version of fasta-formatted files. The header of each sequence should contain an unique protein identifier without whitespace (preceeded by the typical fasta >-symbol):
| ``>Protein_ABC123``
| Between the header and the sequence line can be any number of mutation information lines that start with an <-symbol.
.. note::
  Those lines are not supported by the standard fasta format.
| ``<A23G [Tags-String]``
| The sequence line works identical to the standard fasta format by simply providing an amino acid sequence in one-letter code:
| ``MAGGHKLMRAARAFTP``
| A complete example:
| ``Bla``

.. _structguy_inputs:

Inputs dedicated for further procession with StructGuy
----------------------------------------
Bla
