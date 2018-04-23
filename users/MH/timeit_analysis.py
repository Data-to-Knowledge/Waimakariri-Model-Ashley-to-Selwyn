# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 12/11/2017 2:12 PM
"""

import timeit

number = 1
print('function 1')
print(timeit.timeit('timeit_function()',
                    setup='from users.MH.timeit_function_def import timeit_function',number=number)/number)
print('function 2')
print(timeit.timeit('timeit_function2()',
                    setup='from users.MH.timeit_function_def import timeit_function2',number=number)/number)

