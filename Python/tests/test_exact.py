#!/usr/bin/env python3
import exact_tests
import pytest

def test_fisher():
    pval = exact_tests.fisher([[4, 0], [0, 4]])
    assert pval == pytest.approx(1 / 35)
    pval = exact_tests.fisher([[4, 0], [0, 4]], midp=True)
    assert pval == pytest.approx(1 / 70)
    pval = exact_tests.fisher([[4, 0], [0, 4]], alternative="greater", midp=True)
    assert pval == pytest.approx(1 / 140)
