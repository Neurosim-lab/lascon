: dm/dt = (minf - m)/tau
: input event adds w to m
: when m = 1, or event makes m >= 1 cell fires
: minf is calculated so that the natural interval between spikes is invl

NEURON {
  ARTIFICIAL_CELL INVLF
  RANGE tau, refrac, m, invl, id
  RANGE fflag
  GLOBAL prnum
  : m plays the role of voltage
}

PARAMETER {
  tau = 5 (ms)   <1e-9,1e9>
  invl = 10 (ms) <1e-9,1e9>
  refrac = 2 (ms)
  fflag = 1
  prnum = -1
}

ASSIGNED {
  m
  minf
  t0(ms)
  refractory
  id
}

CONSTRUCTOR {
  VERBATIM 
  if (ifarg(1)) { id= *getarg(1); } else { id= 0; }
  ENDVERBATIM
}

INITIAL {
  if (invl<=refrac) { 
    VERBATIM
    hoc_execerror("invlfire ERROR: invl<refrac.", 0);
    ENDVERBATIM
  }
  minf = 1/(1 - exp(-invl/tau))	: spikes is i
  : m = exp(-invl*scop_random()/tau)
  : m = (1 - minf*(1-m))/m
  m = 0
  t0 = t
  refractory = 0 : 0-integrates input, 1-refractory
  net_send(firetime(), 1)
}

FUNCTION M () { : minf is set so that M will reach 1 after invl with no intervening events
  if (refractory == 0) {
    M = minf + (m - minf)*exp(-(t - t0)/tau)
  }else{
    M = 0 
  }
}

NET_RECEIVE (w) {
  m = M()
  t0 = t
  if (flag == 0) { : external event
    if (refractory == 0) { : otherwise, do nothing
      m = m + w  
      if (m > 1) {
        m = 0
        net_event(t)
        refractory = 1
        net_move(t + refrac)
      } else {
        if (id<prnum) { 
          printf("**** MOVE %g t=%g, to %g (m=%g,%g)\n",id,t,t+firetime(),m,firetime()) }
        net_move(t + firetime())
      }
    }
  } else {
    if (refractory == 0) {
      net_event(t)
      refractory = 1
      net_send(refrac,1)
    } else {
      refractory = 0
      m=0
      net_send(invl-refrac, 1)
    }
  }
}

FUNCTION firetime()(ms) { : m < 0 and minf > 1
  firetime = tau*log((minf-m)/(minf - 1))
  : printf("firetime=%g\n", firetime)
}
