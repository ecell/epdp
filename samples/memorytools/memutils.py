
# Modules that come with Python
import gc
import inspect

# Guppy Heapy tool
from guppy import hpy

# Simple script that outputs memory usage of python
import memoryusage

# The module you want to inspect.
import cluster

def putinfile(output_file, data):

    outFile = open(output_file, 'w')
    outFile.write(str(data))
    outFile.close()

def sort_datav2(class_list, object_list):
    """ Sort objects into class-specific array """

    sorted_data = [[]]
    i = 0

    for a_class in class_list:
           
        for obj in object_list:
            try:
                if isinstance(obj, a_class):
                    sorted_data[i].append(str(obj))
            # when not class
            except:               
                pass 

        i = i + 1 
        sorted_data.append([])

    return sorted_data
            
def check_lengths(array_with_arrays):

    outputarray = []

    total = 0
    i = 0
    for arrays in array_with_arrays:
        outputarray.append(len(arrays))
        print str(i) + ": " + str(len(arrays))
        i = i + 1
        total = total + len(arrays)
    
    return outputarray, total
  
def auto(what):
    """ Automatically runs an analysis. 

    WATCH IT when running from python command line interface! 
    Returns loads of data. Best use is then probably:
    >> d, c, n, t = memutils.auto()    
    """

    classes = giveclasses(what)

    sorted_data = sort_datav2(classes, gc.get_objects())
    lengths_list, total = check_lengths(sorted_data)
    
    return sorted_data, classes, lengths_list, total

def giveclasses():
    """ Returns a list of classes it finds with help from Python's garbage collector."""    

    classlist = []

    for obj in gc.get_objects():
        if inspect.isclass(obj):
            classlist.append(obj)
            
    return classlist

def go(N=100):
    """ (Memory-intensive) code you'd like to analyze. """
    for M in range(N): temp = cluster.single_run(15, False)

def gosmart(filename, C=10, N=100):
    """ Same as go(), but now with some memusage and Heapy stats """
    h = hpy()
    h.setrelheap()

    for G in range(C):    
        for M in range(N): 
            temp = cluster.single_run(15, False)
        hpy().heap().stat.dump(filename)
        print "Current memory usage: " + str(memoryusage.memory())
        
    h.pb(filename)

def gomeasure(C=100, N=10):
    """ Same as go(), but with some memusage stats. """
    print "List of memory usage: "

    for G in range(C):    
        for M in range(N): 
            temp = cluster.single_run(15, False)
        print str(memoryusage.memory())

# NOTES: to call a class via a string:
# C.__dict__["x"]








    
