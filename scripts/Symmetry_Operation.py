#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------------------------------------------
# A moudule store the lattice symmetry operation.
#------------------------------------------------------------------
import numpy as np

def translation(direct_list, cycle):
    """Translate all input direct coordiante in the list.

    Args:
        direct_list (list[list[float]]): The inputed direct coordiante list.
        cycle (float): Translating period

    Returns:
        list[list[float]]: Translated direct coordiante list.
    """
    translation_direct_list = [[direct[0]+cycle[0], 
                                direct[1]+cycle[1], 
                                direct[2]+cycle[2]]
                               for direct in direct_list]
    return translation_direct_list