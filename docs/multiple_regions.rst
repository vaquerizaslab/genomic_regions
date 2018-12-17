.. _multiple-regions:

###########################################
Working with collections of genomic regions
###########################################

.. contents::
   :depth: 2


In the last section, we explored how to work with individual :class:`~GenomicRegion`
objects. This section will demonstrate how to work with lists or collections of regions,
either loaded from file in any of the supported genomic file formats (BED, GFF, BigWig,
Tabix, BEDPE, CSV or tab-delimited tables), or constructed programmatically.


This tutorial assumes you have imported the genomic_regions package like this:

.. code:: python

    import genomic_regions as gr


*******************************************
Loading or constructing RegionBased objects
*******************************************

To load genomic regions from a supported file format (see above), simply use the
:code:`load` function:

.. code:: python

    rb = gr.load("/path/to/file.bed")

:code:`rb` is now a :class:`~RegionBased` object, providing access to a number
of useful methods to work with the regions contained in it.


You can also easily construct a :class:`~RegionBased` object from existing
:class:`~GenomicRegion` lists using :class:`~RegionWrapper`:

.. code:: python

    regions = []
    for chromosome in ['chr1', 'chr2']:
        for start in range(1, 10000000, 100000):
            r = gr.GenomicRegion(chromosome, start, start + 100000,
                                 score=(start - 1)/100000)
            regions.append(r)

    rb = gr.RegionWrapper(regions)

Depending on the size of your region list, the last step could take a long time,
as it constructs an internal data representation of the region list that
facilitates fast searches across intervals. :code:`rb` is now a
:class:`~RegionBased` object.


********************************
Working with RegionBased objects
********************************

~~~~~~~~~~~~~~
Region subsets
~~~~~~~~~~~~~~

The central attribute/method of :class:`~RegionBased` objects is :code:`regions`.
When used as a property, it iterates over all regions in the object:

.. code:: python

    for region in rb.regions:
        print(region)

    # chr1:1-100000
    # chr1:100001-200000
    # ...

If supported by the specific :class:`~RegionBased` subclass (works with most file
types, otherwise a :code:`NotImplementedError` will be thrown) you can access ranges
of, or specific regions using the square bracket notation:

.. code:: python

    for region in rb.regions[0:5]:
        print(region)

    # chr1:1-100001
    # chr1:100001-200001
    # chr1:200001-300001
    # chr1:300001-400001
    # chr1:400001-500001

    print(rb.regions[1])  # chr1:100001-200001

However, the real power of :code:`regions` lies in its double-use as a method.
Without arguments, :code:`regions()` behaves exactly as :code:`regions`. By
providing a region as first argument to :code:`regions()`, you can extract
ranges of regions that overlap with the query:

.. code:: python

    for region in rb.regions('chr1:1-300k'):
        print(region)

    # chr1:1-100001
    # chr1:100001-200001
    # chr1:200001-300001

You can also change the chromosome representation on the fly:

.. code:: python

    for region in rb.regions('chr1:1-300k', fix_chromosome=True):
        print(region)

    # 1:1-100001
    # 1:100001-200001
    # 1:200001-300001

If you are interested in all the chromosomes inside the :class:`~RegionBased`
object, simply use the :code:`chromosomes()` method.


~~~~~~~~~~~~~~
Region binning
~~~~~~~~~~~~~~

If your region objects are associated with scores, i.e. each object has a
:code:`score` attribute with a float value, you can make use of the binning
functions in :class:`~RegionBased` to get binned scores in a defined interval.

The two main methods for this purpose are :code:`binned_regions`, which outputs
:class:`~GenomicRegion` objects, and :code:`region_intervals`, which simply
returns tuples of the form :code:`(start, end, score)`. Other than the
return type, the functions behave in identical fashion, so we are going to focus on
:code:`binned_regions`.

Simply provide :code:`binned_regions` with a genomic interval in the form of a string
or a :class:`~GenomicRegion`, specify the number of :code:`bins`, and you will obtain
equal-sized regions dividing the interval with scores equal to the mean of region
scores falling into each bin, weighted by the size of the associated region.

.. code:: python

    br = rb.binned_regions('chr1', bins=3)

    for region in br:
        print(region, region.score)

    # chr1:1-3333334 16.169996634032984
    # chr1:3333335-6666667 49.3366270640767
    # chr1:6666668-10000001 82.50000000000001

Alternatively, you can specify a :code:`bin_size`:

.. code:: python

    br = rb.binned_regions('chr1', bin_size=3000000)

    for region in br:
        print(region, region.score)

    # chr1:1-3000000 14.499990333423547
    # chr1:3000001-6000000 43.99999032267116
    # chr1:6000001-9000000 73.99999032267117

Note that when choosing a :code:`bin_size` directly, partial bins at the end of the
interval will be omitted.

You can control different aspects of the binning with additional parameters. Most
importantly, you can smooth the scores by choosing a :code:`smoothing_window` size
:code:`n`, which will average scores across :code:`n` neighboring bins up- and
downstream of each bin. I.e. a :code:`smoothing_window` of 2 will average across
5 bins: two to the left, the original bin, and two to the right.
If your data contains NaN values, you can replace them with a fixed value using
:code:`nan_replacement`. On the other hand, you can use :code:`zero_to_nan` to
remove scores of 0 from the calculations.