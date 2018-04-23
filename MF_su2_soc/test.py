from import_file import *

Row=20; Col=Row
g_sbd_l = [1.0,2.0,3.0,4.0]
g_tcur_l = [0.0]#[1.0,2.0,3.0,4.0]

tau = 0.00001

mf.main_cluster(Row, Col, tau, g_sbd_l, g_tcur_l)
