#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import numpy as np
import copy
from matplotlib import pyplot as plt

import scraytrace as sc


# ---- Define tests here
class TestScraytrace(unittest.TestCase):

    def setUp(self):
        pass

    def test_CME(self):
        
        imsize = [256, 256]
        modelid = 54
        modparam = [1.1, 0.52, 5., 0.4, 1e6, 0.2, 0, 0.2, 0.2]
        physics = 0
        DisttoSun_Rsun = 215.
        
        rt = sc.scraytrace(imsize=imsize,
                        frontinteg=True,
                        losrange=[215-20, 215+20],
                        losnbp=100,
                        modelid=modelid, 
                        modparam=modparam, 
                        physics=physics, 
                        fovpix=np.deg2rad(5. / imsize[0]),
                        projtypecode=1,
                        obsang=[0.,0,0],
                        obspos=[0., 0., -DisttoSun_Rsun], 
                        neang=np.deg2rad([30,220,0]), 
                        nbthreads=16, 
                        nbchunks=16,
                        phyparam=[0.4])

        rt.raytrace()
        
        
        palette = copy.copy(plt.cm.gray)
        palette.set_under(color='black')
        palette.set_bad(color='black')

        
        
        rt.dispim(cmap=palette, minmax=(1e-15, 1e-9))

        
        self.assertTrue(True)
        

    def tearDown(self):
        pass

if __name__ == "__main__" :
    unittest.main()
    
    