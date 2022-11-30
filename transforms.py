"""
transforms.py

Created by Vahid Azimi and Mark Zaydman

Purpose: Functions for generating transforms used in PC clonality index

This script contains functions and factory functions to generated 
log, z, and pc transforms using a non-mg cohort data model.

Functions:
	log_transform(INPUT: np.array, inverse: bool = False) -> np.array:
		'''Returns log transformed or inverse log transformed version of input np array'''
	make_ztransform(X: np.array)->Callable[[np.array,bool],np.array]:
		'''Returns z transform function for data model X'''
	make_pctransform(X: np.array, return_decomposition=False)->Callable[[np.array,bool],np.array]:
		'''Returns PC transform function for data model X'''

"""

import numpy as np
import scipy.linalg
from typing import Callable

def log_transform(INPUT: np.array, inverse: bool = False) -> np.array:
	"""Returns log transformed or inverse log transformed version of input np array"""
	OUTPUT=np.copy(INPUT)
	if inverse:
		OUTPUT=np.exp(INPUT)
	else:
		OUTPUT=np.log(INPUT)
	return(OUTPUT)

def make_ztransform(X: np.array)->Callable[[np.array,bool],np.array]:
	"""Returns z transform function for data model X"""
	def z_transform(INPUT: np.array,inverse: bool =False) -> np.array:
		"""Returns the z transformed or inverse z transformed version of the input np array"""
		OUTPUT=np.copy(INPUT)
		if inverse:
			for col in range(INPUT.shape[1]):
				OUTPUT[:,col]=INPUT[:,col]*np.std(X[:,col])+np.mean(X[:,col])
		else:
			for col in range(INPUT.shape[1]):
				OUTPUT[:,col]=(INPUT[:,col]-np.mean(X[:,col]))/np.std(X[:,col])
		return(OUTPUT)
	return(z_transform)

def make_pctransform(X: np.array, return_decomposition=False)->Callable[[np.array,bool],np.array]:
	"""Returns PC transform function for data model X"""
	U,S,Vh=scipy.linalg.svd(X,full_matrices=False)
	def pc_transform(INPUT: np.array,inverse=False) -> np.array:
		"""Returns the PC transformed or inverse PC transformed version of the input np array"""
		if inverse:
			OUTPUT=INPUT@Vh
		else:
			OUTPUT=np.copy(INPUT)
			for row in range(INPUT.shape[0]):
				OUTPUT[row,:]=INPUT[row,:]@Vh.T
		return(OUTPUT)
	if return_decomposition:
		return(pc_transform, U, S, Vh)
	else:
		return(pc_transform)
