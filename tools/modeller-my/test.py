import sys
from modeller import *

sys.stdout = open('test.log', 'w')
print(sys.argv[1])