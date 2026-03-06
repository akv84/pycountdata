"""
countdata — Beta-binomial and inverted beta-binomial tests for count data.

Install:  pip install pycountdata
Import:   import countdata  /  from countdata import bb_test, ibb_test

Functions:
    bb_test     Unpaired beta-binomial test (2+ groups)
    ibb_test    Paired inverted beta-binomial test (2 equal groups)
    normalize   Column-total normalization
    fold_change Per-row fold change

Citations:
    Pham TV et al. (2010) Bioinformatics, 26(3):363-369.
    Pham TV, Jimenez CR (2012) Bioinformatics, 28(18):i596-i602.
"""

from countdata._core import bb_test, ibb_test, normalize, fold_change

__version__ = "1.0.0"
__all__ = ["bb_test", "ibb_test", "normalize", "fold_change"]