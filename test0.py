"""Illustrates the asynchrony of multiple processors"""
from neuron import h
print('all')
pc=h.ParallelContext() # here's the splitting
print(pc.id(),' of ',pc.nhost())
if pc.id()==0: print('Message from rank0 (node0)')
h.quit()
