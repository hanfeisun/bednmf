from numpy import rank
import scipy
from scipy import linalg
import nimfa
from numpy.linalg import pinv

def init_info(model):
    print "init ok"
#    print "Initialized basis matrix\n", model.basis()
#    print "Initialized mixture matrix\n", model.coef()

def factorization(V, rank = 4):
    """
    use nmf to factorize V
    :rtype : (1) the projection matrix (2) the feature vector of V
    """
    fctr = nimfa.mf(V, method="nmf", max_iter=30, rank=rank, update="divergence", objective="div", callback_init=init_info, callback=init_info)
    fctr_res = nimfa.mf_run(fctr)
    print "calculate generized inverse"
    projection = pinv(fctr_res.basis().todense())
    print "inverse finished"
    return {'projection': projection,
            'feature' : (projection * V),
            'basis': fctr_res.basis(),
            'coef': fctr_res.coef().todense()}


#    print "run ok"

