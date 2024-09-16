"""
SPDI unit tests
"""
import pytest
import spdi


def test_right_trim_alleles():
    trimmed = spdi.right_trim_alleles("ATTT", "ATT")
    assert trimmed == ("AT", "A")
    trimmed = spdi.right_trim_alleles("AG", "CG")
    assert trimmed == ("A", "C")


def test_left_trim_alleles():
    trimmed = spdi.left_trim_alleles(100, "AT", "AG")
    assert trimmed == (101, "T", "G")
    trimmed = spdi.left_trim_alleles(100, "T", "G")
    assert trimmed == (100, "T", "G")
    trimmed = spdi.left_trim_alleles(100, "GAT", "GATTG")
    assert trimmed == (102, "T", "TTG")
    trimmed = spdi.left_trim_alleles(100, "CATAA", "CAT")
    assert trimmed == (102, "TAA", "T")
    trimmed = spdi.left_trim_alleles(100, "CATAA", "CAT", over_trim=True)
    assert trimmed == (103, "AA", "")
    trimmed = spdi.left_trim_alleles(100, "GAT", "GATTG", over_trim=True)
    assert trimmed == (103, "", "TG")
    trimmed = spdi.left_trim_alleles(100, "AT", "AG", over_trim=True)
    assert trimmed == (101, "T", "G")
