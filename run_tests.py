import os
import unittest
import sys

minimization_dir = f'{os.getcwd()}/nanome_minimization/'
sys.path.append(minimization_dir)

test_directory = 'tests'
suite = unittest.TestLoader().discover(test_directory)

output = unittest.TextTestRunner(verbosity=1).run(suite)

if output.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)
