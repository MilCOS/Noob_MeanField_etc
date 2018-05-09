import numpy as np
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:25:46 2017
@author: Xuhan
"""
"""
Add square_lattice April 21 2018
"""


def sqr_noflux(Row, Col):
    n = Row*Col
    con = np.zeros([2*n, 2*n], dtype=complex)
    for x in range(Row):
        for y in range(Col):
            i = x*Col + y

            xx = x + 1
            yy = y
            if xx == Row:
                xx = 0
            j = xx*Col + yy
            if i != j:
                for f in [0, 1]:
                    con[i*2+f, j*2+f] = -1.0
                    con[j*2+f, i*2+f] = -1.0

            xx = x
            yy = y + 1
            if yy == Col:
                yy = 0
            j = xx*Col + yy
            if i != j:
                for f in [0, 1]:
                    con[i*2+f, j*2+f] = -1.0
                    con[j*2+f, i*2+f] = -1.0

    return con


def mn_sqr(Row, Col):
    n = Row*Col
    pn = np.ones(n, dtype=int)
    for x in range(Row):
        for y in range(Col):
            i = x * Col + y
            pn[i] = int((-1)**int(x+y))

    return pn
            

def slieb_noflux(Ncell,Nsite,L):
    star = False
    marco = False
    t_hop = np.zeros([2*Nsite,2*Nsite],dtype=float)
    for x in range(L):
        for y in range(L):
            x1 = x * L * Ncell + y * Ncell
            t_hop = get_inside_hopping(x1,t_hop)

            if x == L-1 : star = True
            if y == L-1 : marco = True
            
            t_hop, star, marco = (
                    get_nnext_hopping(x1, t_hop, star, marco, L, Ncell)
                    )
            
    return t_hop

def get_inside_hopping(xin,con):
    
    # vertical
    for j in range(4):
        xt = xin + j*4
        for i in range(3):
            if (j==3)&(i==2): 
                continue
            xx = xt
            xt = xx + 1
            
            con[xx*2,xt*2] = -1.0
            con[xt*2,xx*2] = -1.0
            con[xx*2+1,xt*2+1] = -1.0
            con[xt*2+1,xx*2+1] = -1.0
            
            
    # horizonal
    for i in range(4):
        xt = xin + i
        for j in range(3):
            if (i==3)&(j==2):
                continue
            xx = xt 
            xt = xx + 4
            
            con[xx*2,xt*2] = -1.0
            con[xt*2,xx*2] = -1.0
            con[xx*2+1,xt*2+1] = -1.0
            con[xt*2+1,xx*2+1] = -1.0           
            
    return con
    
def get_nnext_hopping(xin,con,star,marco, L, Ncell):
    
    remainder = (xin//Ncell) % L # modulo of L
    y = remainder # next bolck
    x = ( xin//Ncell - remainder) // L
    
    #vertical
    xb = xin + L * Ncell
    if (star):
        # Boundary return to Row_0: x=0
        xb = (0 * L * Ncell + y * Ncell)
    for i in range(3):
        xx = xin + 3 + 4*i
        xbx = xb + 4*i
        
        con[xx*2,xbx*2] = -1.0
        con[xbx*2,xx*2] = -1.0
        con[xx*2+1,xbx*2+1] = -1.0
        con[xbx*2+1,xx*2+1] = -1.0
        
    #horizonal !!!!
    xr = xin + Ncell
    if (marco):
        # Boundary return to Col_0: y=0
        xr = (x * L * Ncell + 0 * Ncell)
    for j in range(3):
        xx = xin + 3*4 + j
        xrx = xr + j
        
        con[xx*2,xrx*2] = -1.0
        con[xrx*2,xx*2] = -1.0
        con[xx*2+1,xrx*2+1] = -1.0
        con[xrx*2+1,xx*2+1] = -1.0
    
    star = False
    marco = False
    
    return con, star, marco
    
def mn_slieb(Ncell,Nsite):
    mn = np.ones(Ncell,dtype=int)
    mn_table = np.empty(Nsite,dtype=int)
    for j in range(4):
        for i in range(4):
            if (j*i==9): continue
            cat = (-1)**(i+j)
            mouse = j * 4 + i
            mn[mouse] = mn[mouse] * cat
            
    for i in range(Nsite//Ncell):
        # cycling as a Bravis cell
        x1 = i*Ncell
        mn_table[x1: x1+Ncell] = mn[:]
    return mn_table

    
#===== Lieb ======
def lieb_noflux(Ncell,Nsite,L):
    t_hop = np.zeros([2*Nsite,2*Nsite],dtype=float)
    for x in range(L):
        for y in range(L):
            x1 = x*L*3 + y*3
            x2 = x1 + 1
            x3 = x1 + 2
            t_hop[2*x1,2*x2] = -1.0
            t_hop[2*x2,2*x1] = -1.0
            t_hop[2*x1+1,2*x2+1] = -1.0
            t_hop[2*x2+1,2*x1+1] = -1.0

            t_hop[2*x1,2*x3] = -1.0
            t_hop[2*x3,2*x1] = -1.0
            t_hop[2*x1+1,2*x3+1] = -1.0
            t_hop[2*x3+1,2*x1+1] = -1.0
            
            xx = x*L*3 + (y+1)*3 #new cell on the right
            if (y==L-1):
                xx = x*L*3
            t_hop[2*x3,2*xx] = -1.0
            t_hop[2*xx,2*x3] = -1.0
            t_hop[2*x3+1,2*xx+1] = -1.0
            t_hop[2*xx+1,2*x3+1] = -1.0
            
            yy = (x+1)*L*3 + y*3 # new cell on the top
            if (x==L-1):
                yy = y*3
            t_hop[2*x2,2*yy] = -1.0
            t_hop[2*yy,2*x2] = -1.0
            t_hop[2*x2+1,2*yy+1] = -1.0
            t_hop[2*yy+1,2*x2+1] = -1.0

    return t_hop
 
def mn_lieb(Ncell,Nsite):
    mn = np.ones(Nsite)
    for i in range(Nsite):
        if (i%3 != 0):
            mn[i] = -1
    return mn


if __name__ == '__main__':
    f = open('T_hop.txt','w')
    Ncell = 1; L = 4
    Nsite = L*L*Ncell
    Ht = sqr_noflux(L,L)
    mn_table = mn_sqr(L,L)
    for i in range(2*Nsite):
        for j in range(2*Nsite):
            if Ht[i,j] != 0:
                f.write('%i, %i, %3f \n'%(i,j,Ht[i,j]))
    print(mn_table)
    w, v = np.linalg.eigh(Ht)
    print(w)
    f.write(str(mn_table.tolist()))
    f.close()

    
    
    
    
    
