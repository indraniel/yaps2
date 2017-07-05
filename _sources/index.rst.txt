.. toctree::
   :maxdepth: 2
   :caption: Contents:

Introduction
============

``yaps2`` is a loose collection of data processing pipelines meant the `CCDG project <https://www.genome.gov/27563570/>`_.  All of these pipelines are currently designed to be run internally at the `McDonnell Genome Institute (MGI) <https://genome.wustl.edu>`_.

There are currently 4 available pipelines:

#. **postvqsr**

   This is a pipeline to further post-process `VCF files <https://en.wikipedia.org/wiki/Variant_Call_Format>`_ that have been generated from `GATK's <https://software.broadinstitute.org/gatk/documentation/article.php?id=39>`_ `Variant Quality Score Recalibration (VQSR) <https://software.broadinstitute.org/gatk/documentation/article.php?id=39>`_ pipeline.

#. **postvqsr38**

   Similiar to the ``postvqsr`` pipeline, but meant to work on VCF files based on `Human Reference Build version 38 <https://www.ncbi.nlm.nih.gov/grc/human>`_.

#. **Mendelian Inheritance Errors (mie)**

   A niche data processing pipeline to examine the relationship between Mendelian Inheritance Errors to the VQSR score for a selected set of genotypes.

#. **Principal Component Analysis (pca)**

   A niche data processing pipeline to perform principal component analysis on a selected set of genotypes.

#. **Build 38 realignment**

   A data pipeline to realign `BAM <https://samtools.github.io/hts-specs/SAMv1.pdf>`_ / `CRAM <https://samtools.github.io/hts-specs/CRAMv3.pdf>`_ files on Build 38 via |speedseq|_.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`search`

.. |speedseq| replace:: ``speedseq``
.. _speedseq: https://github.com/hall-lab/speedseq
