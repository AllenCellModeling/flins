#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Test myosin
"""

import pytest

import numpy as np
import stress_fiber as sf
import stress_fiber.proteins.myosin as myo


SPAN = 10000
TRACTSPACE = sf.space.TractSpace(1, SPAN)
TRACTS = TRACTSPACE.all_tracts
MYOSINS = [myo.Myosin(np.random.randint(SPAN), t) for t in TRACTS]
# HEADS = [head for myosin in MYOSINS for head in myosin.heads]
# Create bound example
MYOSINS.append(myo.Myosin((np.random.rand() * 0.5 + 0.25) * SPAN, TRACTS[0]))
ACTINS = [sf.proteins.Actin(x, TRACTS[0], 10) for x in MYOSINS[-1].boundaries]
GACTINS = [act.pairs[0] for act in ACTINS]
HEADS = MYOSINS[-1].heads
for gact, head in zip(GACTINS, HEADS):
    head.bs.bind(gact.bs)
    head.state = 1

# FIXME LEFT OFF WORKING ON GETTING THIS TEST SET UP WORKING, WAS THEN GOING TO
# TEST SAMPLE SPRING LENGTHS FOR THE CASE OF FULLY BOUND. WAS THEN GOING TO
# START RE_INTRODUCING KINETICS BACK INTO THE MYOSIN MOTOR.

class TestMyosin:
    @pytest.mark.parametrize("myosin", MYOSINS)
    def test___str__(self, myosin):
        assert not str(myosin).startswith("<")

    @pytest.mark.parametrize("myosin", MYOSINS)
    def test_freely_diffuse(self, myosin):
        xi = myosin.x
        dx = myosin._freely_diffuse()
        xf = xi + dx
        is_x = myosin.x == xf
        at_limits = myosin.heads[0].x_tip == pytest.approx(0) or myosin.heads[
            1
        ].x_tip == pytest.approx(SPAN)
        assert is_x or at_limits

    @pytest.mark.parametrize("myosin", MYOSINS)
    def test__sample_spring_lengths(self, myosin):
        # FIXME Only doing the unbound case currently
        sl_i = [sl for sl in myosin._spring_lengths]
        myosin._sample_spring_lengths()
        sl_f = [sl for sl in myosin._spring_lengths]
        if not myosin.fully_bound:
            assert not any([i == f for i, f in zip(sl_i, sl_f)])
        else:
            assert all([i == f for i, f in zip(sl_i, sl_f)])


    # @pytest.mark.parametrize('myosin', MYOSINS)
    # def test_step(self, myosin):
    #    myosin.step()
