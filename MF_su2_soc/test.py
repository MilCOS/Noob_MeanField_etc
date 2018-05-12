from import_file import *

Row=20; Col=Row
g_sbd_l = [0.01+i*0.2 for i in range(8)]
g_tcur_l = [4.0]

tau = 0.0001
fname = 'sbd_tcur_u'
mf.main_cluster(Row, Col, tau, g_sbd_l, g_tcur_l, fname)
