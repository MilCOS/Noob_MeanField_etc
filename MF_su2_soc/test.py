from import_file import *

Row=20; Col=Row
g_sbd_l = 2.0
g_tcur_l = 0.0

delta = [0 + i*0.01 for i in range(10)]

tau = 0.0001

mf.main_doping(Row, Col, tau, g_sbd_l, g_tcur_l, delta)
