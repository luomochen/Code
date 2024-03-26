#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------
# This script uses the idpp method to interpolate path 
# from initial to final state.
#--------------------------------------------------------
import os
import re
import argparse
from ase.io import vasp
from ase.neb import NEB
from ase.optimize import BFGS

def parse_args():