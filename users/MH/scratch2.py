# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 12/11/2017 2:12 PM
"""

import timeit


print('function 1')
print(timeit.timeit('timeit_function()',setup='from users.MH.scratch import timeit_function',number=10)/10)
print('function 2')
print(timeit.timeit('timeit_function2()',setup='from users.MH.scratch import timeit_function2',number=10)/10)


