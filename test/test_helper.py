'''
Created on 19 mag 2021

@author: lbusoni
'''
import os


def testDataRootDir():
    assert os.path.join(os.path.dirname(__file__), 'data')
