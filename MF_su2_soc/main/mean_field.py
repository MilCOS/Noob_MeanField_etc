import main.lattice_con as lc
import numpy as np
import os
import json


def get_lattice(Row, Col):
    global get_sign
    '''get_lattice():
    get hopping matrix 4x4-1 with spin index
    #1   2   3   ... N
    u d u d u d ... u d
    '''
    Hhop = lc.sqr_noflux(Row, Col)

    get_sign = lc.mn_sqr(Row, Col)

    return Hhop


# ========
def para_solver(rho, ein):
    '''para_slover(rho, ein)
    solve para_new by calculating pho=<c_ic_j> in mean field ansas after diagonalizing process
    para_dic: M_ij, dimer state. We choose the site whoes index is 0;
    para: N_ij, d density wave state. We choose the four sites at the left bottom corner.
          Their indexes are 0, 1, Col, Col+1
    ein[0:4,0:3]: [isite/,UC/x/ysite], for example, for site=0, ein[0,1:]=[1,Col] i.e. neighbor sites' index
                  and column 0 is used to store these four sample sites in ein[:,0]=[0,1,Col,Col+1]
    '''
    para_dic_new = np.empty(8, dtype=complex)

    i = 1
    isite = ein[i, 0]
    ysite, xsite = ein[i, 1:]  # direction: y+, x+
    fisite = isite*2
    fxsite = xsite*2
    fysite = ysite*2

    xdirec = -1j*(
        rho[fisite, fxsite] - rho[fisite+1, fxsite+1]
        - rho[fxsite, fisite] + rho[fxsite+1, fisite+1]
        )/4.0 * get_sign[i]  # i+x
    ydirec = +1j*(
        rho[fisite, fysite] - rho[fisite+1, fysite+1]
        - rho[fysite, fisite] + rho[fysite+1, fisite+1]
        )/4.0 * get_sign[i]  # i+y
    para_new = abs(ydirec)
    if abs(xdirec-ydirec) > 0.1:
        print(xdirec, ydirec, 'N is Not Jerenny: ')
        para_new = (xdirec+ydirec)/2

    for r in range(4):
        isite = ein[r, 0]
        ysite, xsite = ein[r, 1:]  # direction: y+, x+
        fisite = isite*2
        fxsite = xsite*2
        fysite = ysite*2
        xdirec = (
            rho[fisite, fxsite] + rho[fxsite, fisite] +
            rho[fisite+1, fxsite+1] + rho[fxsite+1, fisite+1]
            )  # i+x
        ydirec = (
            rho[fisite, fysite] + rho[fysite, fisite] +
            rho[fisite+1, fysite+1] + rho[fysite+1, fisite+1]
            )  # i+y

        para_dic_new[r*2] = xdirec
        para_dic_new[r*2+1] = ydirec
#        if abs(xdirec-ydirec) > 0.1:
#            print(xdirec, ydirec, 'M is Not Jerenny: '+'-'*r)
#            para_new = (xdirec+ydirec)/2

    return para_dic_new, para_new


def diagonal(ein, para_dic, g1, para, g2):
    '''diagonal(para_ic, U, L)
    Ht is hamilotoian's hopping part;
    para_dic is the mean field parameters of singlet bond order
    g1 is the singlet bond-bond interaction energy;
    para is the mean field parameters of triplet current order;
    g2 is the triplet current-current parameters
    Ht and (para_dic, g1, para, g2) together construct the Hamiltonian which we diagonalise later
    '''
    para_l = para_dic * g1
    para_n = np.array([0, para*g2, -1.0*para*g2])
    M0 = np.zeros([Nsite*2, Nsite*2], dtype=complex)
    N0 = np.zeros([Nsite*2, Nsite*2], dtype=complex)
#    Mu = np.diagflat(np.ones([Nsite*2],dtype=complex))
#    print(Mu)
#    mu = (para_dic[0]**2+para_dic[1]**2+para_dic[2]**2+para_dic[3]**2)*Nsite/4
#    Mu *= mu

    for i in range(0, Nsite):  # rescaled into spin-up channel

        x, y = Aidxy(i, Row, Col)  # get x-y coordinate of site i
        if y+1 != Col:
            j_y = x*Col+y+1
        else:
            j_y = x*Col+0
        if x+1 != Row:
            j_x = (x+1)*Col+y
        else:
            j_x = 0*Col+y

        if (i == j_x)or(i == j_y):
            print("WoW")
            continue
        # spin up
        N0[i*2, j_x*2] = -1j*para_n[get_sign[i]]
        N0[j_x*2, i*2] = +1j*para_n[get_sign[i]]
        N0[i*2, j_y*2] = +1j*para_n[get_sign[i]]
        N0[j_y*2, i*2] = -1j*para_n[get_sign[i]]
        # spin down
        N0[i*2+1, j_x*2+1] = +1j*para_n[get_sign[i]]
        N0[j_x*2+1, i*2+1] = -1j*para_n[get_sign[i]]
        N0[i*2+1, j_y*2+1] = -1j*para_n[get_sign[i]]
        N0[j_y*2+1, i*2+1] = +1j*para_n[get_sign[i]]

        # site: 0
        if (x % 2 == 0) & (y % 2 == 0):
            for f in range(2):
                M0[i*2+f, j_x*2+f] = para_l[0]
                M0[j_x*2+f, i*2+f] = para_l[0]
                M0[i*2+f, j_y*2+f] = para_l[1]
                M0[j_y*2+f, i*2+f] = para_l[1]
        # site: 1
        if (x%2==0)&(y%2!=0):
            for f in range(2):
                M0[i*2+f,j_x*2+f] = para_l[2]
                M0[j_x*2+f,i*2+f] = para_l[2]
                M0[i*2+f,j_y*2+f] = para_l[3]
                M0[j_y*2+f,i*2+f] = para_l[3]
        # site: 2
        if (x%2!=0)&(y%2!=0):
            for f in range(2):
                M0[i*2+f,j_x*2+f] = para_l[4]
                M0[j_x*2+f,i*2+f] = para_l[4]
                M0[i*2+f,j_y*2+f] = para_l[5]
                M0[j_y*2+f,i*2+f] = para_l[5]
        # site: 3
        if (x%2!=0)&(y%2==0):
            for f in range(2):
                M0[i*2+f,j_x*2+f] = para_l[6]
                M0[j_x*2+f,i*2+f] = para_l[6]
                M0[i*2+f,j_y*2+f] = para_l[7]
                M0[j_y*2+f,i*2+f] = para_l[7]

#    print(Ht[:8,:8])
    H = Ht - M0 - N0  #+ Mu
#    print(H[:8,:8])
    w, v = np.linalg.eigh(H)  # w[i]: v[:,i]

    nk = np.zeros([2*Nsite, 2*Nsite], dtype=complex)
    for i in range(Nelec):
        nk[i, i] = 1.0

    rho = np.matmul(np.matmul(v, nk), np.conjugate(v.T))
    para_dic_new, para_new = para_solver(rho, ein)
#    print('rho',(rho[:5,:5]))

    return para_dic_new, para_new, rho


def main_cycle(tau, g1, g2, memo, ein):
    '''main_cycle(tau=0.0001, U, L, memo, flag=0):
    tau is the criteria about shifts between para_old and para_new
    This is the iteration process
    '''
    if memo[0] == 1:
        para_dic_old = abs(np.random.randn(len(memo[1])))  # init para from scartch
        para_old = (np.random.randn(1))
        print('para from init:\n', ' M:', para_dic_old, '\n N:', para_old)
    elif memo[0] > 1:
        print('para from last U:', memo[1:])
        para_dic_old = memo[1]  # init M from last iteration
        para_old = memo[2]  # init N from last iteration

    _tau = 100.0  # init error
    count = 0
    while (_tau > tau):
        para_dic_new, para_new, rho = diagonal(ein, para_dic_old, g1, para_old, g2)
        _tau = np.hstack([abs(para_dic_new - para_dic_old), abs(abs(para_new) - abs(para_old))]).max()
        # print(para_dic_new, para_new)
        if g1 == 0:
            _tau = abs(para_new-para_old).max()
        count += 1
        # print(count,':',para_old)
        para_dic_old = para_dic_new
        para_old = para_new
        if count > 100:
            print('iterations cutoff: ', count)
            break
    return para_dic_new, para_new, rho


def main_cluster(row, col, tau=0.0001, g1list=0.1, g2list=0, fname='tcur'):
    '''
    tau: criteria
    g1list: g_sbd list
    g2list: g_tbd list

    '''
    global Ht, Nsite, Nelec, Row, Col

    Row = row
    Col = col

    Nsite = Row*Col
    Nelec = Nsite

    memo = [0, [i for i in range(8)], 1]  # memory

    ein = np.zeros([4, 3], dtype=int)
    ein[:, 0] = [0, 1, Col+1, Col]
    ein[0, 1:] = [1, Col]  # 0: y,x
    ein[1, 1:] = [2, Col+1]  # 1: y,x
    ein[2, 1:] = [Col+1+1, Col*2+1]  # Col+1: y,x
    ein[3, 1:] = [Col+1, Col*2]  # Col: y,x

    M_U = []
    N_U = []
    g1_t = []
    g2_t = []
    for g_sbd in g1list:
        for g_tcur in g2list:
            Ht = get_lattice(Row, Col)
            memo[0] += 1
            M_sbd, N_tcur, rho = main_cycle(tau, g_sbd, g_tcur, memo, ein)
            memo[1] = M_sbd
            memo[2] = N_tcur

            N_U.append(N_tcur.real)
            M_U.append(M_sbd.real.tolist())
            g1_t.append(g_sbd)
            g2_t.append(g_tcur)
#        print('finished M',g_sbd,g_tcur)
#        print('finished N',g_sbd,g_tcur)

#    fig = plt.figure(figsize=(5,5))
#    ax = fig.add_subplot(111)
#    im = ax.imshow(rho.imag,'Reds')
#    plt.colorbar(im)
#    ax.xaxis.tick_top()
#    ax.set_title('rho')
#    #ax.axis["bottom"].major_ticklabels.set_visible(False)
#    ax.set_anchor('S')
#    plt.show()
    tag = fname

    path = os.getcwd()
    filepath = os.path.join(path, 'results/c1/%s.json' % (tag))
    data = [g1_t, M_U, g2_t, N_U]

    with open(filepath, 'w') as f:
        json.dump(data, f)

    f.close()
    print('g1:\n', g1_t)
    print('g2\n:', g2_t)

    return None


def Aidxy(i, row, col):
    y = i % col
    x = int((i-y)/col)

    return x, y


def main_doping(row, col, tau=0.0001, g1=0.1, g2=0, delta_list=[0]):
    '''
    tau: criteria
    g1: g_sbd *fixed
    g2: g_tbd *fixed
    delta_list: doping of electron
    '''
    global Ht, Nsite, Nelec, Row, Col

    Row = row
    Col = col

    Nsite = Row*Col

    memo = [0, [i for i in range(8)], 1]  # memory

    ein = np.zeros([4, 3], dtype=int)
    ein[:, 0] = [0, 1, Col+1, Col]
    ein[0, 1:] = [1, Col]  # 0
    ein[1, 1:] = [2, Col+1]  # 1
    ein[2, 1:] = [Col+1+1, Col*2+1]  # Col+1
    ein[3, 1:] = [Col+1, Col*2]  # Col

    M_U = []
    N_U = []
    for dope in delta_list:
        Nelec = int(Nsite * (1-dope))  # redefine the electron
        Ht = get_lattice(Row, Col)
        memo[0] += 1
        M_sbd, N_tcur, rho = main_cycle(tau, g1, g2, memo, ein)
        memo[1] = M_sbd
        memo[2] = N_tcur

        N_U.append(N_tcur.real)
        M_U.append(M_sbd.real.tolist())

    path = os.getcwd()
    filepath = os.path.join(path, 'results/c2/gsbd_%i_gtcur_%i_delta.json'%(g1,g2))
    data = [[g1, g2], delta_list, M_U, N_U]

    with open(filepath, 'w') as f:

        json.dump(data, f)

    f.close()
    print('g_sbd,g_tcur: \n', g1, g2)
    print('delta:\n', delta_list)

    return None
