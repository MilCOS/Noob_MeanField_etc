fig, ax = plt.subplots(figsize=(6,4))
Ncell = [3,15][1]
L = 15
Nsite = L*L*Ncell
lname = ['lieb','slieb'][1]
get_sign = -1*mf.get_mn(lname,Ncell,Nsite)

def draw(u1,u2):
    pickle_file = '%2.2f~%2.2f%s_mean_field_L%i_U%2.2f.pkl'%(u1,u2-0.05,lname,L,u2-0.05)
    with open(pickle_file,'rb') as f:
        Mu = []
        Nu = []
        data = pickle.load(f)
        #Ulist = np.arange(u1,u2,0.05)
        Ulist = data.keys()
#        print(data)
        for u in Ulist:
            para_dic = data[u]
            M,N = 0.0,0.0
            for i in range(Ncell):
                Mi = 0.5 * (
                    ( 1 + get_sign[i] * para_dic[i] )
                    - ( 1 - get_sign[i] * para_dic[i] ) 
                    )
                Ni = Mi * get_sign[i]
                M += Mi
                N += Ni
            Mu.append(abs(M))
            Nu.append(abs(N))
    return Ulist,Mu,Nu
c = ['','','','red','black']
nodes = [0.05,0.25,8.0]
#nodes = [0.05,1.0,4.0,8.0]
x,y1,y2 = [],[],[]
for i in range(len(nodes)-1):
    u1,u2 = nodes[i],nodes[i+1]
    Utmp,Mtmp,Ntmp = draw(u1,u2)
    x += Utmp
    y1 += Mtmp
    y2 += Ntmp

ax.plot(x,y2,label='N')
ax.plot(x,y1,label='M')
ax.legend()
ax.set_xlabel('U/t')
ax.set_ylabel('M(N)')
ax.set_title('Mean Field Results')
plt.savefig('mean_field_%s.png'%lname)