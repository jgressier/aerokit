# coding: utf-8
import os
import sys

hadespath = os.path.abspath(__file__)
hadespath = hadespath.split('/')[:-1]
hadespath = '/'+os.path.join(*hadespath)

__all__ = []
__version__ = '1.0.0'
__authors__ = u'Arnaud Grebert, Mikko Folkersma, Gonzalo Saez,\
		Jeremie Gressier, Jaime Vaquero, Julien Bodar,\
		Thibault Bridel-Bertomeu, Thibaut Lunet'
__maintainer__ = "Jeremie Gressier"
__email__ = "Jeremie dot Gressier at isae dot fr"
__status__= "Development"

req_version = (2, 7)
cur_version = sys.version_info

if cur_version < req_version:
	raise Exception("Your Python interpreter %s.%s.%s is too old. "
			"Please consider upgrading." % (cur_version[:3]))


