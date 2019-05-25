import cont_fracs
import pytest
from mpmath import mpf as dec
import gen_consts


@pytest.fixture()
def cf():
    return cont_fracs.ContFrac([[1]], [[1]])


@pytest.mark.parametrize('a,b,num_of_iters,target_val,threshold', [
    ([[2]], [[1, -4, 4]], 550, 4/gen_consts.gen_pi_const()+1, dec('8E-4')),
    ([[6]], [[1, -4, 4]], 500, gen_consts.gen_pi_const()+3, dec('1E-8')),
    ([[1, 2]], [[0, 0, 1]], 200, 4/gen_consts.gen_pi_const(), dec('1E-15')),
    ([[2, 1]], [[0, -1]], 400, gen_consts.gen_e_const()/(gen_consts.gen_e_const()-1), dec('1E-15')),
    ([[1]], [[1]], 400, gen_consts.gen_phi_const(), dec('1E-15')),
    ([[2], [0, 12]], [[-1, 6], [-5, 6]], 400, dec('2')**(1/dec('12'))+1, dec('1E-15')),
])
def test_contfrac(cf, a, b, num_of_iters, target_val, threshold):
    cf.reinitialize(a, b)
    cf.gen_iterations(num_of_iters)
    assert cf.compare_result(target_val) < threshold