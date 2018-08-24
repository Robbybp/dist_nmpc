from __future__ import division
from pyomo.environ import *
from pyomo.dae import * 
from dist_mod_r import m 
from pyomo.core.kernel.expr import exp,sqrt 
from numpy.random import normal as nrn 

__author__ = 'Robert Parker'

# m.nominal_zF[1].set_value(0.4)
# m.nominal_zF[2].set_value(0.2)
# m.nominal_zF[3].set_value(0.4)

solver = SolverFactory('ipopt') 
solver.options['tol'] = 1e-6
solver.options['max_iter'] = 10000

for j in m.set_1_NT:
    m.state_act[j,0].set_value(m.x1set[j,1]) 
    m.state_act[m.NT+j,0].set_value(m.x1set[j,2])
    m.state_act[2*m.NT+j,0].set_value(m.M1set[j])
    m.state_act[3*m.NT+j,0].set_value(m.x2set[j,1])
    m.state_act[4*m.NT+j,0].set_value(m.x2set[j,2])
    m.state_act[5*m.NT+j,0].set_value(m.M2set[j]) 

m.normx[0].set_value( sum(m.state_act[s,0]**2 for s in m.state_set)**(0.5) ) 
m.vtr[-1].set_value(99999.0)
m.ltr[-1].set_value(0.0) 

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
# These are the controls and disturbances:
for k in m.fe:
    m.F[k].set_value(m.nominal_F) 
    m.qF[k].set_value(m.nominal_qF)
    for i in m.set_1_3:
        m.zF[k,i].set_value(m.nominal_zF[i])
    m.D1[k].set_value(m.D1ref.value) 
    m.B1[k].set_value(m.B1ref.value)
    m.LT1[k].set_value(m.LT1ref.value)
    m.VB1[k].set_value(m.VB1ref.value)
    m.D2[k].set_value(m.D2ref.value)
    m.B2[k].set_value(m.B2ref.value)
    m.LT2[k].set_value(m.LT2ref.value)
    m.VB2[k].set_value(m.VB2ref.value) 

# initialize FE vars linearly between IC (set) and setpoint (ref)
for k in m.fep:
    for j in m.set_1_NT:
        for p in m.set_1_2:
            m.x1_0[j,p,k].set_value( m.x1set[j,p].value + k/m.nfe*(m.x1ref[j,p] - m.x1set[j,p].value) ) 
            m.x2_0[j,p,k].set_value( m.x2set[j,p].value + k/m.nfe*(m.x2ref[j,p] - m.x2set[j,p].value) ) 
        m.M1_0[j,k].set_value( m.M1set[j].value + k/m.nfe*(m.M1ref[j] - m.M1set[j].value) ) 
        m.M2_0[j,k].set_value( m.M2set[j].value + k/m.nfe*(m.M2ref[j] - m.M2set[j].value) ) 
for j in m.set_1_NT:
    for p in m.set_1_2:
        m.x1_0[j,p,0].fix(m.x1set[j,p].value)
        m.x2_0[j,p,0].fix(m.x2set[j,p].value) 
    m.M1_0[j,0].fix(m.M1set[j].value)
    m.M2_0[j,0].fix(m.M2set[j].value)

# final collocation point initialization:
for i in m.set_1_NT:
    for f in m.fe:
        for j in m.set_1_2:
            m.x1[i,j,f,3].set_value(m.x1_0[i,j,f+1].value)
            m.x2[i,j,f,3].set_value(m.x2_0[i,j,f+1].value)
for i in m.set_1N:
    for f in m.fe:
        m.M1[i,f,3].set_value(m.M1_0[i,f+1].value)
        m.M2[i,f,3].set_value(m.M2_0[i,f+1].value) 

# initialize other collocation points linearly:
for i in m.set_1_NT:
    for f in m.fe:
        for c in m.set_1_2:
            for j in m.set_1_2:
                m.x1[i,j,f,c].set_value( \
                    m.x1_0[i,j,f].value + \
                    c/3*(m.x1[i,j,f,3].value-m.x1_0[i,j,f].value))
                m.x2[i,j,f,c].set_value( \
                    m.x2_0[i,j,f].value + \
                    c/3*(m.x2[i,j,f,3].value-m.x2_0[i,j,f].value))
for i in m.set_1N:
    for f in m.fe:
        for c in m.set_1_2:
            m.M1[i,f,c].set_value( \
                m.M1_0[i,f].value + \
                c/3*(m.M1[i,f,3].value-m.M1_0[i,f].value))
            m.M2[i,f,c].set_value( \
                m.M2_0[i,f].value + \
                c/3*(m.M2[i,f,3].value-m.M2_0[i,f].value))

# initialize y11 etc via x and alpha...
for i in m.set_1_NTm1:
    for j in m.set_1_2:
        for f in m.fe:
            for c in m.cp:
                m.y_1_1[i,j,f,c].set_value(\
                    m.x1[i,j,f,c].value*m.alpha_m[j,j]) 
                m.y_1_2[i,j,f,c].set_value(\
                    (m.x1[i,1,f,c].value*(m.alpha[1]-1) + \
                    m.x1[i,2,f,c].value*(m.alpha[2]-1))+1) 
                m.y1[i,j,f,c].set_value(\
                    m.y_1_1[i,j,f,c].value/m.y_1_2[i,j,f,c].value)
                m.y_2_1[i,j,f,c].set_value(\
                    m.x2[i,j,f,c].value*m.alpha_m[j,j]) 
                m.y_2_2[i,j,f,c].set_value(\
                    (m.x2[i,1,f,c].value*(m.alpha[1]-1) + \
                    m.x2[i,2,f,c].value*(m.alpha[2]-1))+1) 
                m.y2[i,j,f,c].set_value(\
                    m.y_2_1[i,j,f,c].value/m.y_2_2[i,j,f,c].value)
 
# constant vapor flow rate below feed:
for i in m.set_1_NFm1:
    for f in m.fe:
        m.V1[i,f].set_value(m.VB1[f].value) 
        m.V2[i,f].set_value(m.VB2[f].value)

# constant vapor flow rate above feed
for i in m.set_NF_NTm1:
    for f in m.fe:
        m.V1[i,f].set_value(\
                m.VB1[f].value + (1-m.qF[f].value)*m.F[f].value) 
        m.V2[i,f].set_value(m.VB2[f].value) 

# constant liquid flow rate below feed:
for i in m.set_2_NF:
    for f in m.fe:
        m.L1[i,f].set_value(m.LT1[f].value + m.qF[f].value*m.F[f].value )
        m.L2[i,f].set_value(m.LT2[f].value + m.B1[f].value )

# constant vapor flow rate above feed
for i in m.set_NFp1_NT:
    for f in m.fe:
        m.L1[i,f].set_value(m.LT1[f].value)  
        m.L2[i,f].set_value(m.LT2[f].value) 


# Unneccessary
# Initialize liquid flow rates accordindg to weir formula
#for f in m.fe:
#    for c in m.cp:
#        for i in m.set_2_NF:
#            m.L1[i,f,c].set_value(\
#                m.Kbf*( (m.M1[i,f,c].value - m.Muw + \
#                sqrt((m.M1[i,f,c].value-m.Muw)**2+10**(-8)))/2 )**1.5)
#            m.L2[i,f,c].set_value(\
#                m.Kbf*( (m.M2[i,f,c].value - m.Muw + \
#                sqrt((m.M2[i,f,c].value-m.Muw)**2+10**(-8)))/2 )**1.5)
#        for i in m.set_NFp1_NTm1:
#            m.L1[i,f,c].set_value(\
#                m.Kuf*( (m.M1[i,f,c].value - m.Muw + \
#                sqrt((m.M1[i,f,c].value-m.Muw)**2+10**(-8)))/2 )**1.5)
#            m.L2[i,f,c].set_value(\
#                m.Kuf*( (m.M2[i,f,c].value - m.Muw + \
#                sqrt((m.M2[i,f,c].value-m.Muw)**2+10**(-8)))/2 )**1.5) 
#        m.L1[m.NT,f,c].set_value(m.LT1[f].value)
#        m.L2[m.NT,f,c].set_value(m.LT2[f].value) 
#        
## initialize liquid flow rates at finite elements
#for f in m.fe:
#    for i in m.set_2_NF:
#        m.L1_0[i,f].set_value(\
#            m.Kbf*( (m.M1_0[i,f].value - m.Muw + \
#            sqrt((m.M1_0[i,f].value-m.Muw)**2 +10**(-8)))/2)**1.5)
#        m.L2_0[i,f].set_value(\
#            m.Kbf*( (m.M2_0[i,f].value - m.Muw + \
#            sqrt((m.M2_0[i,f].value-m.Muw)**2 +10**(-8)))/2)**1.5)
#    for i in m.set_NFp1_NTm1:
#        m.L1_0[i,f].set_value(\
#            m.Kuf*( (m.M1_0[i,f].value - m.Muw + \
#            sqrt((m.M1_0[i,f].value-m.Muw)**2 +10**(-8)))/2)**1.5)
#        m.L2_0[i,f].set_value(\
#            m.Kuf*( (m.M2_0[i,f].value - m.Muw + \
#            sqrt((m.M2_0[i,f].value-m.Muw)**2 +10**(-8)))/2)**1.5)
#    m.L1_0[m.NT,f].set_value(m.LT1[f].value)
#    m.L2_0[m.NT,f].set_value(m.LT2[f].value)
#for i in m.set_2_NT:
#    m.L1_0[i,m.nfe].set_value(m.L1ref[i])
#    m.L2_0[i,m.nfe].set_value(m.L2ref[i])
#    m.L1_0[i,0].fix(m.L1_0[i,0].value)
#    m.L2_0[i,0].fix(m.L2_0[i,0].value)

# initialize TC's as linear combinations of component tempteratures
for i in m.set_1_NT: 
    for f in m.fe:
        m.TC1[i,f].set_value(m.x1_0[i,1,f].value*353.3 + \
                             m.x1_0[i,2,f].value*383.8 + \
     (1-m.x1_0[i,1,f].value-m.x1_0[i,2,f].value)*411.5)
        m.TC2[i,f].set_value(m.x2_0[i,1,f].value*353.3 + \
                             m.x2_0[i,2,f].value*383.8 + \
     (1-m.x2_0[i,1,f].value-m.x2_0[i,2,f].value)*411.5)
    m.TC1[i,0].fix(m.TC1[i,0].value)
    m.TC2[i,0].fix(m.TC2[i,0].value)
        
# initialize flow rates and component flow rates with mass balances
for f in m.fe:
    for c in m.cp:
        m.M2dot[1,f,c].set_value(m.L2[2,f].value-m.V2[1,f].value-m.B2[f].value) 

        m.M1dot[1,f,c].set_value(m.L1[2,f].value-m.V1[1,f].value-m.B1[f].value) 

        m.M2dot[m.NT,f,c].set_value(m.V2[m.NT-1,f].value-m.LT2[f].value-m.D2[f].value) 

        m.M1dot[m.NT,f,c].set_value(m.V1[m.NT-1,f].value-m.LT1[f].value-m.D1[f].value) 

#        for i in m.set_NFp1_NTm1:
#            m.M1dot[i,f,c].set_value( m.L1[i+1,f,c].value - m.L1[i,f,c].value + m.V1[i-1,f].value - m.V1[i,f].value )

#            m.M2dot[i,f,c].set_value( m.L2[i+1,f,c].value - m.L2[i,f,c].value + m.V2[i-1,f].value - m.V2[i,f].value )

#        for i in m.set_2_NFm1:
#            m.M1dot[i,f,c].set_value( \
#                m.L1[i+1,f,c].value - m.L1[i,f,c].value + \
#                m.V1[i-1,f].value - m.V1[i,f].value ) 
#            m.M2dot[i,f,c].set_value( \
#                m.L2[i+1,f,c].value - m.L2[i,f,c].value + \
#                m.V2[i-1,f].value - m.V2[i,f].value )

#        m.M2dot[m.NF,f,c].set_value( \
#            m.L2[m.NF+1,f,c].value - m.L2[m.NF,f,c].value + \
#            m.V2[m.NF-1,f].value - m.V2[m.NF,f].value + \
#            m.B1[f].value)
#        m.M1dot[m.NF,f,c].set_value( \
#            m.L1[m.NF+1,f,c].value - m.L1[m.NF,f,c].value + \
#            m.V1[m.NF-1,f].value - m.V1[m.NF,f].value + \
#            m.F[f].value) 

# where is the M2 counterpart to the following constraint?
# ans: it is already given above, this is redundant:
#        m.M1dot[NT,f,c].set_value( \
#           m.V1[NT-1,f].value - m.LT1[f].value - m.D1[f].value ) 

        for j in m.set_1_2:
            m.Mx1dot[1,j,f,c].set_value( \
                m.L1[2,f].value*m.x1[2,j,f,c].value - \
                m.V1[1,f].value*m.y1[1,j,f,c].value - \
                m.B1[f].value*m.x1[1,j,f,c].value )
            m.Mx1dot[m.NT,j,f,c].set_value( \
                m.V1[m.NT-1,f].value*m.y1[m.NT-1,j,f,c].value - \
                m.L1[m.NT,f].value*m.x1[m.NT,j,f,c].value - \
                m.D1[f].value*m.x1[m.NT,j,f,c].value ) 
            m.Mx2dot[1,j,f,c].set_value( \
                m.L2[2,f].value*m.x2[2,j,f,c].value - \
                m.V2[1,f].value*m.y2[1,j,f,c].value - \
                m.B2[f].value*m.x2[1,j,f,c].value )
            m.Mx2dot[m.NT,j,f,c].set_value( \
                m.V2[m.NT-1,f].value*m.y2[m.NT-1,j,f,c].value - \
                m.L2[m.NT,f].value*m.x2[m.NT,j,f,c].value - \
                m.D2[f].value*m.x2[m.NT,j,f,c].value )
        
            for i in m.set_2_NFm1:
                m.x2dot[i,j,f,c].set_value( 1/m.M2_0[i,f].value*(\
                    m.L2[i+1,f].value*m.x2[i+1,j,f,c].value - \
                    m.L2[i,f].value*m.x2[i,j,f,c].value + \
                    m.V2[i-1,f].value*m.y2[i-1,j,f,c].value - \
                    m.V2[i,f].value*m.y2[i,j,f,c].value) ) 

                m.x1dot[i,j,f,c].set_value( 1/m.M1_0[i,f].value*(\
                    m.L1[i+1,f].value*m.x1[i+1,j,f,c].value - \
                    m.L1[i,f].value*m.x1[i,j,f,c].value + \
                    m.V1[i-1,f].value*m.y1[i-1,j,f,c].value - \
                    m.V1[i,f].value*m.y1[i,j,f,c].value )) 
                
            for i in m.set_NFp1_NTm1:
                m.x2dot[i,j,f,c].set_value( 1/m.M2_0[i,f].value*(\
                    m.L2[i+1,f].value*m.x2[i+1,j,f,c].value - \
                    m.L2[i,f].value*m.x2[i,j,f,c].value + \
                    m.V2[i-1,f].value*m.y2[i-1,j,f,c].value - \
                    m.V2[i,f].value*m.y2[i,j,f,c].value) ) 

                m.x1dot[i,j,f,c].set_value( 1/m.M2_0[i,f].value*(\
                    m.L1[i+1,f].value*m.x1[i+1,j,f,c].value - \
                    m.L1[i,f].value*m.x1[i,j,f,c].value + \
                    m.V1[i-1,f].value*m.y1[i-1,j,f,c].value - \
                    m.V1[i,f].value*m.y1[i,j,f,c].value )) 
            
            m.x2dot[m.NF,j,f,c].set_value( 1/m.M2_0[m.NF,f].value*(\
                m.L2[m.NF+1,f].value*m.x2[m.NF+1,j,f,c].value -\
                m.L2[m.NF,f].value*m.x2[m.NF,j,f,c].value + \
                m.V2[m.NF-1,f].value*m.y2[m.NF-1,j,f,c].value - \
                m.V2[m.NF,f].value*m.y2[m.NF,j,f,c].value + \
                m.B1[f].value*m.x1[1,j,f,c].value ))

            m.x1dot[m.NF,j,f,c].set_value( 1/m.M1_0[m.NF,f].value*(\
                m.L1[m.NF+1,f].value*m.x1[m.NF+1,j,f,c].value -\
                m.L1[m.NF,f].value*m.x1[m.NF,j,f,c].value + \
                m.V1[m.NF-1,f].value*m.y1[m.NF-1,j,f,c].value - \
                m.V1[m.NF,f].value*m.y1[m.NF,j,f,c].value + \
                m.F[f].value*m.zF[f,j].value ))
            
            for i in m.set_1N:
                m.x1dot[i,j,f,c].set_value( \
                    (m.Mx1dot[i,j,f,c].value-m.x1[i,j,f,c].value\
                    *m.M1dot[i,f,c].value)/m.M1[i,f,c].value )
                m.x2dot[i,j,f,c].set_value( \
                    (m.Mx2dot[i,j,f,c].value-m.x2[i,j,f,c].value\
                    *m.M2dot[i,f,c].value)/m.M2[i,f,c].value )
### End initialization ###
#p = open('plant.txt','w')

#p.write('%2i'%'k' + '%8s'%'VB1')

# get random noises parameters for disturbances:
for k in m.set_0_Km1:
    m.F_rand[k].set_value( nrn(0,m.F_dev.value) )
    for i in m.set_1_2:
        m.zF_rand[k,i].set_value( nrn(0,m.zF_dev.value) )
    m.zF_rand[k,3].set_value( \
            1 - m.zF_rand[k,1].value - m.zF_rand[k,2].value ) 

for k in m.set_0_Km1:
    m.vkmin1.set_value( m.vtr[k-1] )
    m.lkmin1.set_value( m.ltr[k-1] ) 

    for f in m.fe:
        m.F[f].set_value( m.nominal_F )
        m.pV[f].set_value( m.nominal_pV ) 
        m.qF[f].set_value( m.nominal_qF ) 
        for i in m.set_1_3: m.zF[f,i].set_value( m.nominal_zF[i] ) 

    m.F_r.set_value( \
        m.F[0].value + m.noise_switch*m.F_rand[k].value ) 
    for i in m.set_1_3:
        m.zF_r[i].set_value( \
            m.zF[0,i].value + m.noise_switch*m.zF_rand[k,i].value)
    
    print(k) 

    m.write('dist.nl')
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

    with open('r%i.txt'%k, 'w') as f:
        m.display(ostream=f)
    # m.display()
    
    # record metadata
    m.descepsk[k].set_value(m.desceps.value)
    m.termepsk[k].set_value(m.termeps.value)
    # solver status is a dictionary to be able to store
    # pyomo's discriptive termination condition.
    #m.sol_stat[k] = results.solver.termination_condition
    m.step_cost[k].set_value(m.cost[0].__call__())
    m.vtr[k].set_value(m.vk.value)
    m.ltr[k].set_value(m.lk.value)

    # update initial conditions with plant variables at final FE
    # don't get confused by i and j here
    for j in m.set_1_NT:
        for i in m.set_1_2:
            m.x1set[j,i].set_value(m.x1_r[j,i,3].value)
            m.x2set[j,i].set_value(m.x2_r[j,i,3].value)
        m.M1set[j].set_value(m.M1_r[j,3].value)
        m.M2set[j].set_value(m.M2_r[j,3].value)

        m.state_act[j,       k+1].set_value(m.x1set[j,1].value)
        m.state_act[m.NT+j,  k+1].set_value(m.x1set[j,2].value)
        m.state_act[2*m.NT+j,k+1].set_value(m.M1set[j].value)
        m.state_act[3*m.NT+j,k+1].set_value(m.x2set[j,1].value)
        m.state_act[4*m.NT+j,k+1].set_value(m.x2set[j,2].value)
        m.state_act[5*m.NT+j,k+1].set_value(m.M2set[j].value)

    m.normx[k+1].set_value( sum( m.state_act[s,k+1]**2 for s in m.state_set)**(0.5)) 

    for q in m.fem:
        # Those that have been calculated already:
        m.D1[q].set_value( m.D1[q+1].value ) 
        m.B1[q].set_value( m.B1[q+1].value ) 
        m.LT1[q].set_value( m.LT1[q+1].value ) 
        m.VB1[q].set_value( m.VB1[q+1].value ) 
        m.D2[q].set_value( m.D2[q+1].value ) 
        m.B2[q].set_value( m.B2[q+1].value ) 
        m.LT2[q].set_value( m.LT2[q+1].value ) 
        m.VB2[q].set_value( m.VB2[q+1].value ) 

    # optimal steady state values
    m.D1[m.nfe-1].set_value(m.D1ref.value)
    m.B1[m.nfe-1].set_value(m.B1ref.value)
    m.LT1[m.nfe-1].set_value(m.LT1ref.value)
    m.VB1[m.nfe-1].set_value(m.VB1ref.value)
    m.D2[m.nfe-1].set_value(m.D2ref.value)
    m.B2[m.nfe-1].set_value(m.B2ref.value)
    m.LT2[m.nfe-1].set_value(m.LT2ref.value)
    m.VB2[m.nfe-1].set_value(m.VB2ref.value)

    for i in m.set_1_NT:
        for c in m.cp:
            for j in m.set_1_2:
                for q in m.fem:
                    m.x1[i,j,q,c].set_value( m.x1[i,j,q+1,c].value ) 
                    m.x2[i,j,q,c].set_value( m.x2[i,j,q+1,c].value ) 
                m.x1[i,j,m.nfe-1,c].set_value( m.x1ref[i,j] )
                m.x2[i,j,m.nfe-1,c].set_value( m.x2ref[i,j] ) 
                
    for i in m.set_1_NT:
        for j in m.set_1_2:
            for q in m.fe:
                m.x1_0[i,j,q].set_value( m.x1_0[i,j,q+1].value ) 
                m.x2_0[i,j,q].set_value( m.x2_0[i,j,q+1].value ) 
            m.x1_0[i,j,m.nfe].set_value( m.x1ref[i,j] )
            m.x2_0[i,j,m.nfe].set_value( m.x2ref[i,j] )
            m.x1_0[i,j,0].unfix()
            m.x2_0[i,j,0].unfix() 
       
    for i in m.set_1N:
        for c in m.cp:
            for q in m.fem:
                m.M1[i,q,c].set_value( m.M1[i,q+1,c].value )
                m.M2[i,q,c].set_value( m.M2[i,q+1,c].value ) 
            m.M1[i,m.nfe-1,c].set_value( m.M1ref[i] )
            m.M2[i,m.nfe-1,c].set_value( m.M2ref[i] )


    for i in m.set_1_NT:
        for q in m.fe:
            m.M1_0[i,q].set_value( m.M1_0[i,q+1].value )
            m.M2_0[i,q].set_value( m.M2_0[i,q+1].value )
        m.M1_0[i,m.nfe].set_value( m.M1ref[i] )
        m.M2_0[i,m.nfe].set_value( m.M2ref[i] )
        m.M1_0[i,0].unfix()
        m.M2_0[i,0].unfix()

    for i in m.set_1_NTm1:
        for c in m.cp:
            for j in m.set_1_2:
                for q in m.fem: 
                    m.y_1_1[i,j,q,c].set_value( m.y_1_1[i,j,q+1,c].value ) 
                    m.y_1_2[i,j,q,c].set_value( m.y_1_2[i,j,q+1,c].value ) 
                    m.y_2_1[i,j,q,c].set_value( m.y_2_1[i,j,q+1,c].value ) 
                    m.y_2_2[i,j,q,c].set_value( m.y_2_2[i,j,q+1,c].value ) 
                    # Why is y1 initialized like this?
                    m.y1[i,j,q,c].set_value( m.y_2_2[i,j,q+1,c].value ) 
                    m.y2[i,j,q,c].set_value( m.y2[i,j,q+1,c].value )
            # and y_1_1 etc., for that matter ... 
            m.y_1_1[i,j,m.nfe-1,c].set_value( m.y1ref[i,j] ) 
            m.y_1_2[i,j,m.nfe-1,c].set_value( m.y1ref[i,j] ) 
            m.y_2_1[i,j,m.nfe-1,c].set_value( m.y2ref[i,j] ) 
            m.y_2_2[i,j,m.nfe-1,c].set_value( m.y2ref[i,j] ) 
            m.y1[i,j,m.nfe-1,c].set_value( m.y1ref[i,j] )
            m.y2[i,j,m.nfe-1,c].set_value( m.y2ref[i,j] )

    for i in m.set_1_NTm1:
        for q in m.set_1_nfem2:
           m.V1[i,q].set_value( m.V1[i,q+1].value )
           m.V2[i,q].set_value( m.V2[i,q+1].value )
        m.V1[i,m.nfe-1].set_value( m.V1ref[i] )
        m.V2[i,m.nfe-1].set_value( m.V2ref[i] )
        
    for i in m.set_2_NT:
#        for c in m.cp:
        for q in m.set_1_nfem2:
            m.L1[i,q].set_value( m.L1[i,q+1].value )
            m.L2[i,q].set_value( m.L2[i,q+1].value )
        m.L1[i,m.nfe-1].set_value( m.L1ref[i] ) 
        m.L2[i,m.nfe-1].set_value( m.L2ref[i] ) 

#    for i in m.set_2_NT:
#        for q in m.set_1_nfem2: 
#            m.L1_0[i,q].set_value( m.L1_0[i,q+1].value )     
#            m.L2_0[i,q].set_value( m.L2_0[i,q+1].value )     
#        m.L1_0[i,m.nfe].set_value( m.L1ref[i] )
#        m.L2_0[i,m.nfe].set_value( m.L2ref[i] )
#        m.L1_0[i,0].unfix()
#        m.L2_0[i,0].unfix()

    for i in m.set_1_NT:
        for q in m.fe:
            m.TC1[i,q].set_value( m.TC1[i,q+1].value )
            m.TC2[i,q].set_value( m.TC2[i,q+1].value )
        m.TC1[i,m.nfe].set_value( m.TC1ref[i] ) 
        m.TC2[i,m.nfe].set_value( m.TC2ref[i] )
        m.TC1[i,0].unfix()
        m.TC2[i,0].unfix()


    

m.avg_cost.set_value( sum( m.step_cost[i].value for i in m.set_0_Km1 )/m.K)

    
