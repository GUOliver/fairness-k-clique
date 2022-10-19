import random
import sys
'''

this file used to generator attributes for unlimited node

notice:
the commond line is:
python3 attribute_gen.py number_of_node number_of_attributes

'''

if(len(sys.argv[1])==0):
    print("plz input number of node\n")
    exit

if(len(sys.argv[2])==0):
    print("plz input number of attributes\n")
    exit
N_node = int(sys.argv[1])
N_attribute = int(sys.argv[2])

f = open("attribute_gen.txt","w")
for i in range(0,N_node):
    f.write(str(i)+' '+str(random.randint(0,N_attribute-1))+'\n')

f.close()


