from __future__ import division
from pyomo.environ import *
from pyomo.dae import * 
from ss_mod import m 
from pyomo.core.kernel.expr import exp,sqrt 
from numpy.random import normal as nrn 
from pyutilib.services import TempfileManager
import csv

__author__ = 'Robert Parker'

# m.nominal_zF[1].set_value(0.4)
# m.nominal_zF[2].set_value(0.2)
# m.nominal_zF[3].set_value(0.4)

TempfileManager.tempdir="/home/rparker1/dist_nmpc/steady_state"
solver = SolverFactory('ipopt') 
solver.options['tol'] = 1e-6
solver.options['max_iter'] = 10000

#M1set = {}
#M2set = {}
#M1init = {}
#M2init = {}
#x1init = {}
#x2init = {}
#for i in m.set_1_NT:
#    M1set[i] = m.M1set[i].value
#    M2set[i] = m.M2set[i].value
#for i in m.set_1_NT:
#    m.M1set[i].set_value( m.M1ref[i].value + 0.1*(M1set[i] - m.M1ref[i].value) )
#    m.M2set[i].set_value( m.M2ref[i].value + 0.1*(M2set[i] - m.M2ref[i].value) )
#m.M1set[1].set_value( 3.0 )
#m.M2set[1].set_value( 3.0 )
#m.M1ref[1].set_value( 3.0 )
#m.M2ref[1].set_value( 3.0 )
#m.M1set[m.NT].set_value( 3.0 )
#m.M2set[m.NT].set_value( 3.0 )
#m.M1ref[m.NT].set_value( 3.0 )
#m.M2ref[m.NT].set_value( 3.0 )

# keep track of initial conditions even after *set is replaced:
#for i in m.set_1_NT:
#    M1init = m.M1set[i].value
#    M2init = m.M2set[i].value
#    for j in m.set_1_2:
#        x1init[(i,j)] = m.x1set[i,j].value
#        x2init[(i,j)] = m.x2set[i,j].value

m.vtr.set_value(99999.0)
m.ltr.set_value(0.0) 

# (maybe unnecessary) initialization
m.vk.set_value(0.0)
m.lk.set_value(0.0)

##### 
# Here, skip the part of the .run file defining weights, 
# at least until I start doing eNMPC
#####

############################
### Initialize first NLP ###
############################

# when to use .value and when not to appears 
# inconsistent: 
m.F.set_value(m.nominal_F) 
m.qF.set_value(m.nominal_qF)
for i in m.set_1_3:
    m.zF[i].set_value(m.nominal_zF[i])
m.D1.set_value(m.D1ref.value) 
m.B1.set_value(m.B1ref.value)
m.LT1.set_value(m.LT1ref.value)
m.VB1.set_value(m.VB1ref.value)
m.D2.set_value(m.D2ref.value)
m.B2.set_value(m.B2ref.value)
m.LT2.set_value(m.LT2ref.value)
m.VB2.set_value(m.VB2ref.value) 

# initialize FE vars linearly between IC (set) and setpoint (ref)
# ^unnecessary because FE's no longer exist 
# "_0" vars have no meaning

# "differential state" initialization:
for i in m.set_1_NT:
    for j in m.set_1_2:
        m.x1[i,j].set_value(m.x1ref[i,j])
        m.x2[i,j].set_value(m.x2ref[i,j])
    m.M1[i].set_value(m.M1ref[i].value)
    m.M2[i].set_value(m.M2ref[i].value) 

# initialize y11 etc via x and alpha...
for i in m.set_1_NTm1:
    for j in m.set_1_2:
       m.y_1_1[i,j].set_value(\
           m.x1[i,j].value*m.alpha_m[j,j]) 
       m.y_1_2[i,j].set_value(\
           (m.x1[i,1].value*(m.alpha[1]-1) + \
           m.x1[i,2].value*(m.alpha[2]-1))+1) 
       m.y1[i,j].set_value(\
           m.y_1_1[i,j].value/m.y_1_2[i,j].value)
       m.y_2_1[i,j].set_value(\
           m.x2[i,j].value*m.alpha_m[j,j]) 
       m.y_2_2[i,j].set_value(\
           (m.x2[i,1].value*(m.alpha[1]-1) + \
           m.x2[i,2].value*(m.alpha[2]-1))+1) 
       m.y2[i,j].set_value(\
           m.y_2_1[i,j].value/m.y_2_2[i,j].value)
 
# constant vapor flow rate below feed:
for i in m.set_1_NFm1:
    m.V1[i].set_value(m.VB1.value) 
    m.V2[i].set_value(m.VB2.value)

# constant vapor flow rate above feed
for i in m.set_NF_NTm1:
    m.V1[i].set_value(\
            m.VB1.value + (1-m.qF.value)*m.F.value) 
    m.V2[i].set_value(m.VB2.value) 

# Initialize liquid flow rates accordindg to weir formula
for i in m.set_2_NF:
    m.L1[i].set_value(\
        m.Kbf*( (m.M1[i].value - m.Muw + \
        sqrt((m.M1[i].value-m.Muw)**2+10**(-8)))/2 )**1.5)
    m.L2[i].set_value(\
        m.Kbf*( (m.M2[i].value - m.Muw + \
        sqrt((m.M2[i].value-m.Muw)**2+10**(-8)))/2 )**1.5)
for i in m.set_NFp1_NTm1:
    m.L1[i].set_value(\
        m.Kuf*( (m.M1[i].value - m.Muw + \
        sqrt((m.M1[i].value-m.Muw)**2+10**(-8)))/2 )**1.5)
    m.L2[i].set_value(\
        m.Kuf*( (m.M2[i].value - m.Muw + \
        sqrt((m.M2[i].value-m.Muw)**2+10**(-8)))/2 )**1.5) 
m.L1[m.NT].set_value(m.LT1.value)
m.L2[m.NT].set_value(m.LT2.value) 
        
# initialize TC's as linear combinations of component tempteratures
for i in m.set_1_NT: 
    m.TC1[i].set_value(m.x1[i,1].value*353.3 + \
                         m.x1[i,2].value*383.8 + \
 (1-m.x1[i,1].value-m.x1[i,2].value)*411.5)
    m.TC2[i].set_value(m.x2[i,1].value*353.3 + \
                         m.x2[i,2].value*383.8 + \
 (1-m.x2[i,1].value-m.x2[i,2].value)*411.5)
        
# initialize flow rates and component flow rates with mass balances

### End initialization ###
#p = open('plant.txt','w')

#err_x = {}
#err_M = {}
#
#
#writer = csv.writer(open('data.csv','a',newline=''), delimiter=',')
#writer.writerow(['k','cp','err_x','err_M','F','M1_F','M2_F','M1_NTm1','M2_NTm1','LT1','VB1','D1','B1','LT2','VB2','D2','B2'])

# get random noises parameters for disturbances:
# ^delete anything referencing noise 

m.F.set_value( m.nominal_F )
m.pV.set_value( m.nominal_pV ) 
m.qF.set_value( m.nominal_qF ) 
for i in m.set_1_3: m.zF[i].set_value( m.nominal_zF[i] ) 

#m.write('dist.row')
#m.write('dist.col')
# ^does not work
#m.write('dist_ss.nl',io_options={'symbolic_solver_labels':True})

# m.const1_r.display()
# m.y_1_1_r.display()
# m.x1_r.display()
# m.L1.display()
# m.V1.display()
# m.M1dot.display()
# m.const11.display()
# m.const8.display()
# m.const12.display()
# for var in m.component_map(Var).itervalues():
#     var.display()
# for con in m.component_map(Constraint).itervalues():
#     con.display() 
# m.display()
results = solver.solve(m,tee=True) 

m.write('dist.nl',io_options={'symbolic_solver_labels':True})

with open('m_ss.txt', 'w') as f:
    m.display(ostream=f)
# m.display()

# record metadata
# solver status is a dictionary to be able to store
# pyomo's discriptive termination condition.
#m.sol_stat[k] = results.solver.termination_condition
m.step_cost.set_value(m.cost.__call__())


    
