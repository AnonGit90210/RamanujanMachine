import lhs_evaluators
import gen_consts
# from decimal import Decimal as dec
from mpmath import mpf as dec, mp
from utils import MathOperations
import pytest

import decimal
# decimal.getcontext().prec=40
mp.dps = 40

@pytest.mark.parametrize('target_constant', [
    gen_consts.gen_pi_const(),
    gen_consts.gen_e_const(),
    gen_consts.gen_phi_const(),
])
@pytest.mark.parametrize('u, l, c, d', [
    (1, 2, 3, 4),
    (10, 12, 59, 9),
    (0, 1, 0, 1),
    (-9, 18, 7, -3)
])
def test_ulcd_evaluator(u, l, c, d, target_constant):
    ulcd_evaluator = lhs_evaluators.ULCDEvaluator((u, l, c, d), target_constant)
    expected_value = (u/target_constant+target_constant/l+c)/d
    assert abs(ulcd_evaluator.get_val() - expected_value) < dec('1E-38')
    ulcd_evaluator.add_int(3)
    expected_value += 3
    assert abs(ulcd_evaluator.get_val() - expected_value) < dec('1E-38')
    ulcd_evaluator.flip_sign()
    expected_value *= -1
    assert abs(ulcd_evaluator.get_val() - expected_value) < dec('1E-38')


@pytest.mark.parametrize('target_constant', [
    gen_consts.gen_pi_const(),
    gen_consts.gen_e_const(),
    gen_consts.gen_phi_const(),
])
@pytest.mark.parametrize('numerator_coeffs, denominator_coeffs', [
    ([0, 1, 1, 5], [5]),
    ([1, -8, 17], [2, 3, 4]),
    ([1, -8, 17], [19, -12, 5, 20]),
    ([1], [19, 1])
])
def test_rationalfunc_evaluator(numerator_coeffs, denominator_coeffs, target_constant):
    rationalfunc_evaluator = lhs_evaluators.RationalFuncEvaluator((numerator_coeffs, denominator_coeffs, 0),
                                                                  target_constant)
    expected_value = MathOperations.subs_in_polynom(numerator_coeffs, target_constant) / \
                     MathOperations.subs_in_polynom(denominator_coeffs, target_constant)
    assert abs(rationalfunc_evaluator.get_val() - expected_value) < dec('1E-38')
    rationalfunc_evaluator.add_int(3)
    expected_value += 3
    assert abs(rationalfunc_evaluator.get_val() - expected_value) < dec('1E-38')
    rationalfunc_evaluator.flip_sign()
    expected_value *= -1
    assert abs(rationalfunc_evaluator.get_val() - expected_value) < dec('1E-38')


@pytest.mark.parametrize('u_range, l_range, c_range, d_range, results', [
    ([-1, 2], [-1, 2], [-1, 2], [-1, 2], {
        (-1, -1, 0, -1),
        (-1, -1, 0, 1),
        (-1, 1, 0, -1),
        (-1, 1, 0, 1),
        (0, -1, 0, -1),
        (0, -1, 0, 1),
        (0, 1, 0, -1),
        (0, 1, 0, 1),
        (1, -1, 0, -1),
        (1, -1, 0, 1),
        (1, 1, 0, -1),
        (1, 1, 0, 1),
    }),
    ([1, 2], [1, 2], [-3, 5], [-2, 2], {
        (1, 1, -1, -2),
        (1, 1, 0, -2),
        (1, 1, 0, -1),
        (1, 1, 0, 1),
        (1, 1, 1, -2),
    })
])
def test_ulcd_enumerator(u_range, l_range, c_range, d_range, results):
    ulcd_enum = lhs_evaluators.ULCDEnumerator([u_range, l_range, c_range, d_range], target_constant=1)
    ulcds_generator = ulcd_enum.generator()
    assert { ulcd_eval.get_params() for ulcd_eval in ulcds_generator }.symmetric_difference(results) == set()


@pytest.mark.parametrize('numerator_coeffs, denominator_coeffs, force_numerator_greater, target_constant, results', [
    ([[-1, 2], [-1, 2]], [[-1, 2], [-1, 2]], False, dec(5), {
        ((-1, -1), (-1, -1)),
        ((-1, -1), (-1, 0)),
        ((-1, -1), (-1, 1)),
        ((-1, -1), (0, -1)),
        ((-1, -1), (0, 1)),
        ((-1, -1), (1, -1)),
        ((-1, -1), (1, 0)),
        ((-1, -1), (1, 1)),
        ((-1, 0), (-1, -1)),
        ((-1, 0), (-1, 1)),
        ((-1, 0), (0, -1)),
        ((-1, 0), (0, 1)),
        ((-1, 0), (1, -1)),
        ((-1, 0), (1, 1)),
        ((-1, 1), (-1, -1)),
        ((-1, 1), (-1, 0)),
        ((-1, 1), (-1, 1)),
        ((-1, 1), (0, -1)),
        ((-1, 1), (0, 1)),
        ((-1, 1), (1, -1)),
        ((-1, 1), (1, 0)),
        ((-1, 1), (1, 1)),
        ((0, -1), (-1, -1)),
        ((0, -1), (-1, 0)),
        ((0, -1), (-1, 1)),
        ((0, -1), (0, -1)),
        ((0, -1), (0, 1)),
        ((0, -1), (1, -1)),
        ((0, -1), (1, 0)),
        ((0, -1), (1, 1)),
        ((0, 1), (-1, -1)),
        ((0, 1), (-1, 0)),
        ((0, 1), (-1, 1)),
        ((0, 1), (0, -1)),
        ((0, 1), (0, 1)),
        ((0, 1), (1, -1)),
        ((0, 1), (1, 0)),
        ((0, 1), (1, 1)),
        ((1, -1), (-1, -1)),
        ((1, -1), (-1, 0)),
        ((1, -1), (-1, 1)),
        ((1, -1), (0, -1)),
        ((1, -1), (0, 1)),
        ((1, -1), (1, -1)),
        ((1, -1), (1, 0)),
        ((1, -1), (1, 1)),
        ((1, 0), (-1, -1)),
        ((1, 0), (-1, 1)),
        ((1, 0), (0, -1)),
        ((1, 0), (0, 1)),
        ((1, 0), (1, -1)),
        ((1, 0), (1, 1)),
        ((1, 1), (-1, -1)),
        ((1, 1), (-1, 0)),
        ((1, 1), (-1, 1)),
        ((1, 1), (0, -1)),
        ((1, 1), (0, 1)),
        ((1, 1), (1, -1)),
        ((1, 1), (1, 0)),
        ((1, 1), (1, 1)),
    }),
    ([[-1, 2], [-1, 2]], [[-1, 2], [-1, 2]], False, dec(1), {
        ((-1, -1), (-1, -1)),
        ((-1, -1), (-1, 0)),
        ((-1, -1), (0, -1)),
        ((-1, -1), (0, 1)),
        ((-1, -1), (1, 0)),
        ((-1, -1), (1, 1)),
        ((-1, 0), (-1, -1)),
        ((-1, 0), (0, -1)),
        ((-1, 0), (0, 1)),
        ((-1, 0), (1, 1)),
        ((0, -1), (-1, -1)),
        ((0, -1), (-1, 0)),
        ((0, -1), (0, -1)),
        ((0, -1), (0, 1)),
        ((0, -1), (1, 0)),
        ((0, -1), (1, 1)),
        ((0, 1), (-1, -1)),
        ((0, 1), (-1, 0)),
        ((0, 1), (0, -1)),
        ((0, 1), (0, 1)),
        ((0, 1), (1, 0)),
        ((0, 1), (1, 1)),
        ((1, 0), (-1, -1)),
        ((1, 0), (0, -1)),
        ((1, 0), (0, 1)),
        ((1, 0), (1, 1)),
        ((1, 1), (-1, -1)),
        ((1, 1), (-1, 0)),
        ((1, 1), (0, -1)),
        ((1, 1), (0, 1)),
        ((1, 1), (1, 0)),
        ((1, 1), (1, 1))
    }),
    ([[-1, 2], [-1, 2]], [[-1, 2], [-1, 2]], True, dec(1), {
        ((-1, -1), (-1, -1)),
        ((-1, -1), (-1, 0)),
        ((-1, -1), (0, -1)),
        ((-1, -1), (0, 1)),
        ((-1, -1), (1, 0)),
        ((-1, -1), (1, 1)),
        ((0, -1), (-1, -1)),
        ((0, -1), (-1, 0)),
        ((0, -1), (0, -1)),
        ((0, -1), (0, 1)),
        ((0, -1), (1, 0)),
        ((0, -1), (1, 1)),
        ((0, 1), (-1, -1)),
        ((0, 1), (-1, 0)),
        ((0, 1), (0, -1)),
        ((0, 1), (0, 1)),
        ((0, 1), (1, 0)),
        ((0, 1), (1, 1)),
        ((1, 1), (-1, -1)),
        ((1, 1), (-1, 0)),
        ((1, 1), (0, -1)),
        ((1, 1), (0, 1)),
        ((1, 1), (1, 0)),
        ((1, 1), (1, 1))
    })
])
def test_rationalfunc_enumerator(numerator_coeffs, denominator_coeffs, force_numerator_greater, target_constant,
                                 results):
    rationalfunc_enum = lhs_evaluators.RationalFuncEnumerator([numerator_coeffs, denominator_coeffs,
                                                               force_numerator_greater], target_constant)
    rationalfunc_generator = rationalfunc_enum.generator()
    assert { rationalfunc_eval.get_params()[:2]
             for rationalfunc_eval in rationalfunc_generator }.symmetric_difference(results) == set()


def test_rationalfunc_canonalize():
    rationalfunc_eval = lhs_evaluators.RationalFuncEvaluator([[0, 3*1, 3*2, 1], [0, 1, 2], 0], dec(5))
    rationalfunc_eval.canonalize_params()
    rationalfunc_eval.get_params() == ([0, 0, 1], [1, 2], 3)