from jitcdde import y

# instantiate jitcdde state variables
c_f1 = lambda tau=None: y(0, tau) if tau is not None else y(0)
c_f2 = lambda tau=None: y(1, tau) if tau is not None else y(1)
c_t1 = lambda tau=None: y(2, tau) if tau is not None else y(2)
c_c1 = lambda tau=None: y(3, tau) if tau is not None else y(3)
c_c2 = lambda tau=None: y(4, tau) if tau is not None else y(4)
c_m = lambda tau=None: y(5, tau) if tau is not None else y(5)

n = lambda tau=None: y(6, tau) if tau is not None else y(6)
C1 = lambda tau=None: y(7, tau) if tau is not None else y(7)
rho = lambda tau=None: y(8, tau) if tau is not None else y(8)
