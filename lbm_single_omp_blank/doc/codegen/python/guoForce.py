#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:31:12 2019

@author: diogo
"""

import numpy as np
import sympy as sp

# Weights and vector for the d3q19 lattice
W0 = sp.Rational(1,3)
W1 = sp.Rational(1,18)
W2 = sp.Rational(1,36)

sp.init_printing()

i   = sp.Idx("i", (0,18) )
alfa = sp.Idx("\\alpha", (0,2) )
beta = sp.Idx("\\beta", (0,2) )

w = sp.IndexedBase("w",i)
v = sp.IndexedBase("v",alfa)
F = sp.IndexedBase("F",beta)
c = sp.IndexedBase("c", (alfa,i) )

cs = sp.Symbol("c_s")
tau = sp.Symbol("tau")
lbd = sp.Symbol("lambda")

# Fguo = w[i] * ( 1- 1/(2*tau)) *( (c[alfa,i]  - v[alfa]) / (cs**2) )  * c[alfa,i]
Fguo = w[i] * ( 1- 1/(2*tau))  *( (c[alfa,i]  - v[alfa]) / (cs**2)  +  c[beta,i] * v[beta] * c[alfa,i] / (cs**4) )*F[alfa]
FguoOdd =  w[i] * ( 1- 1/(2*tau))  *(  - sp.Sum( v[alfa]*F[alfa], alfa) / (cs**2)  + sp.Sum( c[beta,i] * v[beta] , beta) * sp.Sum( c[alfa,i]*F[alfa] , alfa) / (cs**4) )
FguoEven = w[i] * ( 1- 1/(2*tau))  *( (c[alfa,i]*F[alfa] ) / (cs**2)  )


d3q19_c =  np.array( [ [  0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 ] , \
             [  0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0, -1,  1, -1,  1 ] , \
             [  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1 ]  ]  )

d3q19_w =  [ sp.Rational(1,3) for i in range(0,1) ] + \
           [ sp.Rational(1,18) for i in range(1,7) ] + \
           [ sp.Rational(1,36) for i in range(7,19) ]
           
d3q19_cs = 1/sp.sqrt(3)


for i in range(1,18,2):
    print( np.dot( np.array(  [ sp.Symbol("Fx") , sp.Symbol("Fy"), sp.Symbol("Fz") ] ) , d3q19_c[:,i] ) *  np.dot( np.array(  [ sp.Symbol("vx")  , sp.Symbol("vy") , sp.Symbol("vz") ] ) , d3q19_c[:,i] ) )

  
#
#vm = np.array(  sp.symbols("vmx vmy vmz") )
#Fm = np.array(  sp.symbols("Fmx Fmy Fmz") )
#Fmv = sp.Symbol("Fmv")
#
#w  = np.array( [ W0, W1, W1, W1, W1, W1, W1, W2, W2, W2, W2, W2, W2, W2, W2, W2, W2, W2, W2 ] )
#
#for i in range(1,len(cx),2):
#    c = np.array( [cx[i], cy[i] , cz[i] ] )
#    print( "even = "  + sp.printing.ccode( w[i]* ( np.dot( Fm, c)  *  np.dot( vm, c) - Fmv ) ) + ";" )
#    print( "odd = " + sp.printing.ccode( w[i]* np.dot( vm, c) )  + ";")
#    print( "F[" + str(i) + "] = even + odd;" )
#    print( "F[" + str(i+1) + "] = even - odd" + ";\n" )