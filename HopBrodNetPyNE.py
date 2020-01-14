# based on netpyne/doc/source/code/tut7.py
from netpyne import sim, specs
cfg, par = specs.SimConfig(), specs.NetParams()

par.popParams['hop'] = {'cellType': 'PYR', 'cellModel': 'HH', 'numCells': 50}
par.cellParams['PYR'] = {'conds': {'cellModel': 'HH', 'cellType': 'PYR'},  
                          'secs': {'soma' : {'geom': {'diam': 18.8, 'L': 18.8, 'Ra': 123.0}, 'mechs': {'hh':{}}, 'vinit':-70}}}
par.synMechParams['exc'] = {'mod': 'Exp2Syn', 'tau1': 0.1, 'tau2': 1.0, 'e': 0}
par.synMechParams['inh'] = {'mod': 'Exp2Syn', 'tau1': 0.1, 'tau2': 1.0, 'e': -80}
par.stimSourceParams['bkg'] = {'type': 'NetStim', 'rate': 50, 'noise': 0.5}
par.stimTargetParams['bkg->all'] = {'source': 'bkg', 'conds': {'pop': 'hop'}, 'weight': 0.1, 'delay': 1, 'synMech': 'exc'}
par.connParams['hop->hop'] = { 'preConds': {'pop': 'hop'}, 'postConds': {'pop': 'hop'}, 'weight': 0.2, 'synMech': 'inh', 'delay': 5}

cfg.duration, cfg.dt = 500, 0.025
cfg.analysis['plotRaster']={'syncLines': True}; # cfg.analysis['plotTraces']={'include': [1]}; cfg.analysis['plot2Dnet']=True
sim.createSimulateAnalyze(netParams = par, simConfig = cfg)
