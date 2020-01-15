from neuron import h
import sys, pylab as plt, glob, numpy as np, pickle as pkl
from itertools import permutations as permut
h.load_file('stdrun.hoc')
plt.ion()

useHHCell = True # toggle for using HH cells vs INVLF ArtCells
tvec, ind, = h.Vector(), h.Vector()
cells, nclist, ncell = [], [], 100

class Cell:
    def __init__ (self, _id, _x=0, _y=0, _z=0):
        self.pp = h.INVLF()
        self.pp.tau, self.pp.invl = 10, 20
        self.gid, self.x, self.y, self.z = _id, _x, _y, _z

class HHCell:
    def __init__ (self, _id, _x=0, _y=0, _z=0, record=True):
        self.soma = h.Section(name='soma', cell=self)
        self.soma.insert('hh')
        self.stim = h.IClamp(0.5, sec=self.soma)
        self.stim.amp, self.stim.delay, self.stim.dur = 10, 0, 1e9
        self.syn = h.Exp2Syn(0.5, sec=self.soma)
        self.syn.tau1, self.syn.tau2, self.syn.e = 0.1, 1, -80
        self.gid, self.x, self.y, self.z = _id, _x, _y, _z
        if record: 
          self.vvec = h.Vector()
          self.vvec.record(self.soma(0.5)._ref_v)

def createnet (ncell=100):
    global ncells, cells, nclist # need global since want to set these
    ncells = ncell
    cells = [HHCell(i) if useHHCell else Cell(i,i,0,0) for i in range(ncell)]
    nclist = [h.NetCon(x[0].soma(0.5)._ref_v, x[1].syn, sec=x[0].soma) if useHHCell else h.NetCon(x[0].pp,x[1].pp) for x in permut(cells,2)] 
    for ce in cells:
        nc= h.NetCon(ce.soma(0.5)._ref_v, None, sec=ce.soma) if useHHCell else h.NetCon(ce.pp, None)
        nc.record(tvec, ind, ce.gid)
        nclist.append(nc)
      
firststep = True
def setparams (w=None, delay=None, tau=None, low=None, high=None): 
    """params for both netcons and cells"""
    global firststep
    for n in nclist: 
        if w is not None: n.weight[0] = w if useHHCell else -w # neg wts are needed for ArtCell (neg driving force for HHCell)
        if delay is not None: n.delay=delay
    if firststep or (low and high): # only set this up using defaults the 1st time
        firststep = False 
        if low is not None: low, high = 10, 11
        for ce,ran,start in zip(cells, np.random.uniform(low, high, ncells), np.random.uniform(0, 30, ncells)): 
            if useHHCell:
                ce.stim.amp = ran
                ce.stim.delay = start
            else: 
                ce.pp.invl = ran
                if tau is not None: ce.tau = tau

# run options -- run() autorun() mpirun()
def run (ncell=100, tstop=500):
  h.tstop = tstop
  createnet(ncell)
  setparams(w=0.0, delay=2, tau=10, low=10, high=12)
  h.run()

def mpirun (ncell=100, tstop=500, w=0.4):
    '''run under mpi; need to set final mpirun() command in this file so it runs!
       eg `mpiexec -np 4 --oversubscribe nrniv -python -mpi net.py`    '''   
    global pc, idhost, nhost, ncells, cells, gidlist, nclist, useHHCell, data, gather
    useHHCell, ncells, h.tstop, = True, ncell, tstop
    pc = h.ParallelContext() # nodes split in different directions
    pc.set_maxstep(10)
    idhost, nhost = pc.id(), pc.nhost()
    print(idhost, " of ",nhost)
    gidlist = range(idhost, ncells, nhost) # round robin
    cells = [HHCell(i,record=True) for i in gidlist]
    for ce in cells:
        pc.set_gid2node(ce.gid, idhost)
        nc = h.NetCon(ce.soma(0.5)._ref_v, None, -20, 0, 0, sec=ce.soma) # presynaptic monitor for connections
        pc.cell(ce.gid, nc, 1)
        pc.spike_record(ce.gid, tvec, ind) # Record spikes of this cell
        del nc # don't need it since isn't providing an actual connection
    nclist = [pc.gid_connect(pre, post.syn) for pre in range(ncells) for post in cells if pre!=post.gid] # pre's are all cells, posts are on this node
    setparams(w=w, delay=2, tau=10, low=10, high=12)
    h.finitialize()
    pc.psolve(h.tstop)  # actually run the sim in all nodes
    data = [None]*nhost
    data[0] = {'tvec': tvec, 'ind': ind}
    pc.barrier()
    gather=pc.py_alltoall(data)
    pc.barrier()

outdict = {}
def autorun (ls = [0.5,0.3,0.1,0.01,0.001]):
    global outdict
    createnet(100)
    setparams(w=0.0, delay=2, tau=10, low=10, high=12)
    for w in ls:
        setparams(w=w)
        h.run()
        s = syncer()
        outdict[w]={'ind': np.array(ind).astype(int), 'tvec': np.array(tvec), 'sync': s}
        print(w, round(s,2), end=" ")

# plotting and analysis
gl, g = [], None
def plotiv (raster=True, cl=[]):
    global gl,g
    if raster:
        if not g: g = h.Graph()
        g.erase_all()
        ind.mark(g, tvec,"O",4,2,1)
        for t in tvec: 
            g.beginline(3,1)
            g.line(t,0)
            g.line(t,ncells)
        g.exec_menu("View = plot")
    if cl:
        gl = [h.Graph() for i in range(len(cl))]
        for c,g in zip(cl,gl): 
            cells[c].vvec.plot(g)
            g.exec_menu("View = plot")
        
def plotpy (tvec=tvec, ind=ind, synclines=True):
    """same as plotiv() but somehow much slower"""
    plt.scatter(tvec,ind)    
    if synclines:
        for spkt in np.array(tvec): plt.plot((spkt, spkt), (0, ncells), 'r-', linewidth=0.1)

def syncer (vec=tvec,width=1):
    '''measures how well spikes "fill up" the time; assumes spike times in tvec, tstop; param: width
    syncer doesn't take account of prob of overlaps due to too many spikes stuffed into too little time'''
    t0=-1; cnt=0
    for tt in vec:
        if tt>=t0+width: 
            t0, cnt = tt, cnt+1
    return 1.0-cnt/(h.tstop/width)

# sample usage
'''
useHHCell = True and False # optionally use HH cells instead of default INVLF cells
run()
plotiv()
setparams(w=0.04)
h.run()
plotiv()
'''

'''
autorun()
plt.plot(outdict.keys(), [x['sync'] for x in outdict.values()])
plotpy(outdict[-0.1]['tvec'], outdict[-0.1]['ind'], False)
'''

mpirun(w=0.4)
# plt.scatter(gather[2]['tvec'],gather[2]['ind'])
