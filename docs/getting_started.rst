###############
Getting started
###############

.. contents::
   :depth: 2


********
Overview
********

The genomic_regions package aims to simplify working with genomic region / interval data
by providing a common interface that lets you access a wide selection of file types and
formats for handling genomic region data - all using the same syntax.

For example, simply do

.. code:: python

    import genomic_regions as gr

    f = gr.load('/path/to/file.<bed|gff|bw|...>')
    for region in f.regions('chr1:10Mb-25Mb'):
        print(region)

    # chr1:10000000-11000000
    # chr1:12000000-12000000
    # ...

to print all regions on chromosome 1 between 10 and 25 megabases. This works with BED,
GFF, BigWig, Tabix, and other file formats, and you can easily extend it yourself!


************
Installation
************

The simplest way to install genomic_regions is via pip:

.. code:: bash

   pip install genomic_regions

and that should be all you need! If you are not the owner of the Python installation,
try:

.. code:: bash

   pip install --user genomic_regions

You can also directly download the genomic_regions source code from Github by cloning its repository.
The installation is then done via setup.py:

.. code:: bash

   git clone http://www.github.com/vaquerizaslab/genomic_regions
   cd genomic_regions
   pip install .

genomic_regions can now be imported as a Python (2.6+, 3.4+) module
