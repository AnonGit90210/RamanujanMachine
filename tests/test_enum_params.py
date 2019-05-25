import enum_params
from gen_consts import gen_pi_const
from postprocfuncs import EVALUATED_POSTPROC_FUNCS, POSTPROC_FUNCS, INVERSE_POSTPROC_PAIRS, POSTPROC_FUNCS_LATEX
from enum_poly_params import BasicEnumPolyParams
from lhs_evaluators import LHS_TYPES
import pytest

# TODO: finish implementing tests. Only some methods are currently tested.

@pytest.fixture
def mitm():
    # postproc_funcs_filter: id, 1/x, sqrt(x)
    return enum_params.MITM(target_generator=gen_pi_const, target_name='pi', postproc_funcs_filter=(0, 1,  10),
                            trunc_integer=True, hashtable_prec=15, ab_poly_class=BasicEnumPolyParams,
                            ab_poly_special_params=None, enum_only_exp_conv=True, num_of_iterations=300, prec=50,
                            postproc_funcs=EVALUATED_POSTPROC_FUNCS, postproc_funcs_text=POSTPROC_FUNCS,
                            postproc_funcs_text_inverse=INVERSE_POSTPROC_PAIRS,
                            postproc_funcs_latex=POSTPROC_FUNCS_LATEX)


@pytest.mark.parametrize('range_a, range_b, results', [
    ([[[1, 2], [2, 3]]], [[[0, 1], [0, 1], [1, 2]]], {
        ((((((1, 2),), ((0, 0, 1),)), 0)),): '0.273',
        ((((((1, 2),), ((0, 0, 1),)), 1)),): '0.7853',
        ((((((1, 2),), ((0, 0, 1),)), 10)),): '0.128',
    })
])
def test_gen_hashtable(mitm, range_a, range_b, results):
    mitm.build_hashtable(range_a, range_b)
    assert set(results.keys()).symmetric_difference(set([ tuple(v) for v in mitm.dec_hashtable.values() ])) == set()
    for k, values in mitm.dec_hashtable.items():
        print('\n\n', k, '\n\n')
        for (a, b), postproc_func_id in values:
            assert k.startswith(results[(((a, b), postproc_func_id),)])



@pytest.mark.parametrize('range_a, range_b, lhs_classes, lhs_enum_params, result_clicks', [
    ([[[1, 2], [2, 3]]], [[[0, 1], [0, 1], [1, 2]]], LHS_TYPES['ulcd'], [[0, 1], [4, 5], [0, 1], [1, 2]], {
        ((((1, 2),), ((0, 0, 1),)), (0, 4, 0, 1), 1)
    })
])
def test_find_clicks(mitm, range_a, range_b, lhs_classes, lhs_enum_params, result_clicks):
    mitm.build_hashtable(range_a, range_b)
    mitm.find_clicks(lhs_classes, lhs_enum_params)
    filtered_clicks = [(ab, lhs_eval.get_params(), post_func_ind)
                       for ab, lhs_eval, post_func_ind, _ in mitm.get_filtered_params()]
    assert set(filtered_clicks).symmetric_difference(result_clicks) == set()