from jitcdde import y

# generator function, so variables can be defined in arbitrary order below
def counter():
    count = 0
    while True:
        yield count
        count += 1

# instantiate generator object 
c = counter()

# instantiate jitcdde state variables 
c_f1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
c_f2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
c_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
c_c1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
c_c2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
c_m = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_fh1_f1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh1_f2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh1_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh1_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh1_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_fh2_f1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh2_f2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh2_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh2_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_fh2_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_ch1_c1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch1_c2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch1_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch1_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch1_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_ch2_c1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch2_c2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch2_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch2_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_ch2_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_hwf1_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf1_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf1_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf1_w1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf1_w2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_hwf2_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf2_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf2_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf2_w1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwf2_w2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_hwc1_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc1_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc1_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc1_w1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc1_w2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

hx_hwc2_h1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc2_h2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc2_t1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc2_w1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
hx_hwc2_w2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

n = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
C1 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
C2 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
C3 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
C4 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
C5 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
C6 = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))
rho = lambda tau=None: y(next(c), tau) if tau is not None else y(next(c))

