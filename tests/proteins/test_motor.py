#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Test the simple motor
"""

import pytest

import numpy as np
import stress_fiber as sf
from stress_fiber.proteins.motor import Motor


def motor_tract():
    span = 10000
    tractspace = sf.space.TractSpace(0, span)
    t = tractspace.all_tracts[0]
    t.rand = lambda: span * (np.random.rand() * 0.5 + 0.25)
    return t


def unbound_motor():
    t = motor_tract()
    act = sf.proteins.Actin(t.rand(), t, length=100)
    return Motor(t.rand(), t)


def left_side_bound_motor():
    t = motor_tract()
    motor = Motor(t.rand(), t)
    act = sf.proteins.Actin(motor.x, t, length=100)
    motor.heads[0].bs.bind(act.pairs[0].bs)
    motor.heads[0].state = 1
    return motor


def right_side_bound_motor():
    t = motor_tract()
    motor = Motor(t.rand(), t)
    act = sf.proteins.Actin(motor.locs[1], t, length=100)
    motor.heads[1].bs.bind(act.pairs[0].bs)
    motor.heads[1].state = 1
    return motor


def both_sides_weak_motor():
    t = motor_tract()
    motor = Motor(t.rand(), t)
    acts = [sf.proteins.Actin(x, t, length=100) for x in motor.locs]
    gacts = acts[0].pairs[0], acts[1].pairs[0]
    motor.heads[0].bs.bind(gacts[0].bs)
    motor.heads[0].state = 1
    motor.heads[1].bs.bind(gacts[1].bs)
    motor.heads[1].state = 1
    return motor


def both_sides_strong_motor():
    motor = both_sides_weak_motor()
    for head in motor.heads:
        head.state = 2
    return motor


motor_list = [
    unbound_motor(),
    left_side_bound_motor(),
    right_side_bound_motor(),
    both_sides_weak_motor(),
    both_sides_strong_motor(),
]

head_list = [head for motor in motor_list for head in motor.heads]


class TestMotor:
    @pytest.mark.parametrize("motor", motor_list)
    def test___str__(self, motor):
        assert not str(motor).startswith("<")

    @pytest.mark.parametrize("motor", motor_list)
    def test_step(self, motor):
        motor.locs
    
    @pytest.mark.parametrize("motor", motor_list)
    def test_step(self, motor):
        motor.step()

    @pytest.mark.parametrize("motor", motor_list)
    def test_head_access(self, motor):
        for i in (0, 1):
            assert motor.heads[i].other_head is motor.heads[i ^ 1]

    @pytest.mark.parametrize("motor", motor_list)
    def test_force(self, motor):
        pass

    @pytest.mark.parametrize("motor", motor_list)
    def test_energy(self, motor):
        pass

    @pytest.mark.parametrize("motor", motor_list)
    def test_boundaries(self, motor):
        pass
