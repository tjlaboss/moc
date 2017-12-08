# Flux mapping
#
# A module containing a mutable dctionary type for fluxes

import collections

class FixedDict(collections.MutableMapping):
	"""A wrapper for standard dictionaries and mappings.
		
	This class allows for mutable flux dicts
	Copied from @bereal:
		https://stackoverflow.com/a/14816620
	
	Parameters:
	-----------
	data:           original dictionary to wrap
	"""
	def __init__(self, data):
		self.__data = data
	
	def __str__(self):
		return str(self.__data)
	
	def __len__(self):
		return len(self.__data)
	
	def __iter__(self):
		return iter(self.__data)
	
	def __setitem__(self, k, v):
		#if k not in self.__data:
		#	raise KeyError(k)
		self.__data[k] = v
	
	def __delitem__(self, k):
		raise NotImplementedError
	
	def __getitem__(self, k):
		return self.__data[k]
	
	def __contains__(self, k):
		return k in self.__data

