#######################################
Working with individual genomic regions
#######################################

.. contents::
   :depth: 2


Genomic intervals, or regions, as we call them from here on, are represented
by a :class:`~api/regions/GenomicRegion` object. This object has attributes
commonly associated with genomic regions, such as "chromosome", "start" and "end",
but can in principle have arbitrary attributes, including scores, labels, and other
useful properties describing the region. There is no restriction regarding the types
of attributes - any valid Python object can be associated with a
:class:`~genomic_regions.GenomicRegion`.

This tutorial assumes you have imported the genomic_regions package like this:

.. code:: python

    import genomic_regions as gr


************************
Creating genomic regions
************************

You can create a genomic region by calling the :class:`~GenomicRegion` constructor:

.. code:: python

    region = gr.GenomicRegion(chromosome='chr1', start=1, end=1000)

    # or to simplify
    region = gr.GenomicRegion('chr1', 1, 1000)

:code:`start` and :code:`end` must be of type int.
The :code:`strand` attribute also has special restrictions. It can either be a str
(:code:`"+"`, :code:`"-"`, :code:`"."`), an int (:code:`-1`, :code:`+1`) or :code:`None`.

.. code:: python

    region = gr.GenomicRegion('chr1', 1, 1000, strand='+')

Other attributes have no restrictions, but we recommend that :code:`score` be a float,
to show the expected behavior when working with the region later on.


.. code:: python

    region = gr.GenomicRegion('chr1', 1, 1000, strand='+',
                              score=1.765, name="my region",
                              my_list=[1, 2, 3, 4])

You can also add attributes to the region later on by using the :code:`set_attribute`
method:

.. code:: python

    region = gr.GenomicRegion('chr1', 1, 1000, strand='+')
    region.set_attribute("my_dict", {'a': 1, 'b': 2})

We advise the use of :code:`set_attribute` rather than the builtin :code:`setattr`
:code:`region.my_dict = {'a': 1, 'b': 2}`, as some processing is done to the
key, value pair by :class:`~GenomicRegion` for compatibility.

Finally, you can quickly create :class:`~GenomicRegion` objects from strings using the
:code:`as_region` convenience function:

.. code:: python

    region = gr.as_region('chr1:1-1000:+')

The region string should have the format :code:`<chromosome>:<start>-<end>[:<strand>]`.
:code:`start` and :code:`end` can use common abbreviations for kilo- and megabases,
support decimal and thousand separators, and are case-insensitive,
so writing :code:`gr.as_region('chr12:12500000-18000000')` is the same as
:code:`gr.as_region('chr12:12.5Mb-18Mb')` and :code:`gr.as_region('chr12:12.5mb-18mb')` and
:code:`gr.as_region('chr12:12,500,000-18,000,000')`.


**********************
Genomic region methods
**********************

~~~~~
Basic
~~~~~

The :class:`~GenomicRegion` object comes loaded with useful attributes and methods,
most of which are self-explanatory:

.. code:: python

    len(region)  # returns the size of the region in base pairs
    region.center  # returns the base (or fraction of base) at the center of the region
    region.five_prime  # returns the starting base at the 5' end of the region
    region.three_prime  # returns the starting base at the 3' end of the region
    region.is_forward()  # True if strand is '+' or '+1'
    region.is_reverse()  # True if strand is '-' or '-1'
    region.attributes  # return all attribute names in this region object
    region.copy()  # return a shallow copy of this region
    region.to_string()  # return a region identifier string describing the region

The :code:`strand` attribute returns an integer (or :code:`None`, if no strand is set).
To obtain a string, use the method :code:`strand_string`, which returns one of
:code:`+`, :code:`-`, or :code:`.`.

~~~~~~~~~~~~~~~~~~~~
Modifying the region
~~~~~~~~~~~~~~~~~~~~

Some methods are provided that modify the underlying region.

:code:`region.expand` changes the size of the region on the chromosome, either by an
absolute amount in base pairs (using any of the parameters :code:`absolute`,
:code:`absolute_left`, or :code:`absolute_right`), or relative, as a fraction of the
current region size (:code:`relative`, :code:`relative_left`, or :code:`relative_right`).
By default, these actions return a modified copy of the original region, but you can
modify the region in place using :code:`copy=True`.

.. code:: python

    region = gr.as_region('chr12:12.5Mb-18Mb')
    print(region)  # chr12:12500000-18000000
    new_region = region.expand(absolute='1mb')
    print(new_region)  # chr12:11500000-19000000
    print(region)  # chr12:12500000-18000000
    region.expand(relative=1.5, copy=False)
    print(region)  # chr12:4250000-26250000


You can also easily move a region on the same chromosome by adding or subtracting base
pairs.

.. code:: python

    region = gr.as_region('chr12:12.5Mb-18Mb')
    new_region = region + 1000000
    print(new_region)  # chr12:13500000-19000000

Some databases store chromosome names with the 'chr' prefix, others without. You can use
the method :code:`fix_chromosome` to switch between chromosome formats:


.. code:: python

    region = gr.as_region('chr12:12.5Mb-18Mb')
    new_region = region.fix_chromosome()
    print(new_region)  # 12:12500000-18000000


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Relationship to other regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can easily check if a region overlaps with another region:

.. code:: python

    region = gr.as_region('chr12:12.5Mb-18Mb')
    region.overlaps('chr12:11Mb-13Mb')  # True
    region.overlaps('chr12:11Mb-11.5Mb')  # False
    region.overlaps('chr1:11Mb-13Mb')  # False

Similarly, you can get the extent of the overlap in base pairs:

.. code:: python

    region = gr.as_region('chr12:12.5Mb-18Mb')
    region.overlap('chr12:11Mb-13Mb')  # 500000
    region.overlap('chr12:11Mb-11.5Mb')  # 0


**********
Next steps
**********

Next, we will see how to work with lists and collections of regions in
:ref:`multiple-regions`.
