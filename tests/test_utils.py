import utils
import pytest
import time


def test_measure_runtime_measures():
    mt = utils.MeasureRuntime()
    time_measured = []
    sleeping_times = [0.3, 0.2, 0.5]
    mt.start_measure()
    for sleeping_time in sleeping_times:
        time.sleep(sleeping_time)
        time_measured.append(mt.measure_time())
    assert mt.get_total_runtime() - sum(sleeping_times) < 0.1
    assert time_measured[-1] == mt.get_last_measured_time()
    assert all([abs(sleeping_t - measured_t) < 0.03 for sleeping_t, measured_t in zip(sleeping_times, time_measured)])
    assert mt.is_started()


@pytest.fixture
def math_op():
    return utils.MathOperations()


@pytest.mark.parametrize('m,n,res', [
    (1480, 276, 4),
    (3, 0, 3),
    (1, 5, 1),
    (6, -4, -2)
])
def test_math_gcd(math_op, m, n, res):
    assert math_op.gcd(m, n) == res


@pytest.mark.parametrize('poly_coeffs,x,res', [
    ([1, 2, 3, 4, 5], 7, 13539),
    ([1, 2, 3, 4, 5], -5, 2691),
    ([1, -2, 0, -4, 5], 3, 292),
    ([1, -2, 0, -4, 5], -5, 3636)
])
def test_math_substitute_in_poly(math_op, poly_coeffs, x, res):
    assert math_op.subs_in_polynom(poly_coeffs, x) == res