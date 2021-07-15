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
        """Raytrace and display a GCS CM model"""
        
        imsize = [256, 256]
        modelid = 54
        # model parameters: see models51to60.cpp, CModel54::dumpDefaultParamForIDL
        # 
        # ("rb","2.55","dist to bottom of structure","Rsun"));    
        # ("alpha","0.52","angle between Oy and leg axis","rad"));    
        # ("h","6.","leg height","Rsun"));    
        # ("kappa","0.4","aspect ratio","")); 
        # ("nemin","1e6","Ne","cm^-3"));  
        # ("thick","0.1","Thickness of the skeleton axis.","Rsun"));
        # ("neaxis","0.","Electron density in the skeleton or axis","cm^-3"));
        # ("stiffness","0.","NOT USED",""));  
        # ("skinsigmain","0.1","inner sigma",""));    
        # ("skinsigmafr","0.1","front sigma",""));
        alpha = 0.52
        h = 5.
        kappa = 0.4
        modparam = [1.1, alpha, h, kappa, 1e6, 0.2, 0, 0.2, 0.2]
        physics = 0
        DisttoSun_Rsun = 215.

        # compute CME leading edge height
        leadingEdgeHeight = sc.model54_calcLeadingEdgeHeight(h, kappa, alpha)
        
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
                        phyparam=[0.58])

        rt.raytrace()
        
        palette = copy.copy(plt.cm.gray)
        palette.set_under(color='black')
        palette.set_bad(color='black')

        rt.dispim(cmap=palette, minmax=(1e-15, 1e-9))
    
        print("Leading edge height {0} Rsun".format(leadingEdgeHeight))

    
        self.assertTrue(True)
        




    def tearDown(self):
        pass

if __name__ == "__main__" :
    unittest.main()
    
    