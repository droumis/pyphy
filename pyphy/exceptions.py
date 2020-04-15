# -*- coding: utf-8 -*-

"""
pyphy.exceptions
~~~~~~~~~~~~~~~~~~~
This module contains the set of Pyphy's custom exceptions.
"""
    
class Error(Exception):
    def __init__(self, *args, **kwargs):
        """Base class for exceptions in pyphy"""
        super(Error, self).__init__(*args, **kwargs)

class SourceError(Error):
    def __init__(self, *args, source='', **kwargs):
        print(f'An error occurred with source: {source}')
        super(SourceError, self).__init__(*args, **kwargs)

        
class MissingKwargsError(Error):
    def __init__(self, *args, include_kwargs=[], **kwargs):
        print(f'missing kwargs: {include_kwargs}')
        super(MissingKwargsError, self).__init__(*args, **kwargs)