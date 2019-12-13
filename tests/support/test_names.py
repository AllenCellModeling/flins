#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test names.
"""

import flins.support.names as names


sample_names = [names.unique_name() for i in range(20)]


def test_unique_name():
    """Check uniqueness"""
    assert len(sample_names) == len(list(set(sample_names)))


def test_name_to_uuid():
    """Check that we can roundtrip our names"""
    for name in sample_names:
        roundtrip = names.unique_name(names.name_to_uuid(name))
        assert name == roundtrip
