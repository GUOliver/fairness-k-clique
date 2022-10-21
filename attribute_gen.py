import random
import sys
'''

this file used to generator attributes for unlimited node

notice:
the commond line is:
python3 attribute_gen.py number_of_node number_of_attributes

'''

if (len(sys.argv)<3):
    print("[usage] python3 attribute_gen.py number_of_node number_of_attributes\n")
    exit(0)

N_node = int(sys.argv[1])
N_attribute = int(sys.argv[2])

f = open("attribute.txt", "w")
for i in range(0, N_node):
    f.write(str(i) + ' ' + str(random.randint(0, N_attribute - 1)) + '\n')

f.close()


