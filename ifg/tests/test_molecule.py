"""Pytest file for ifg module testing"""

from chem.molecule import Molecule

def test_mol():
    """Place a test molecule to view if the output works"""
    assert Molecule("c1ccc2c(c1)cccc2C#Cc1ccccc1C#Cc1cccc2ccccc12")

def test_fg():
    """Place a test functional group to view if the output works"""
    assert Molecule("[R]C(=O)[R]", name="Ketone", type="fg")
