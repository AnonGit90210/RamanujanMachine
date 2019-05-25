import pytest
import enum_poly_params


@pytest.fixture
def basic_poly():
    return enum_poly_params.BasicEnumPolyParams(num_of_iterations=10, enum_only_exp_conv=False,
                                                avoid_int_roots=True, should_gen_contfrac=True,
                                                avoid_zero_b=True, threshold=None, prec=80)


@pytest.mark.parametrize('poly,should_invert', [
    ([0, 0, 10, -5, 0, 0], 0),
    ([10, -5], 0),
    ([-1, 3, 1], 1)
])
def test_basic_does_have_integer_roots(basic_poly, poly, should_invert):
    assert basic_poly._does_have_int_positive_roots(poly) ^ should_invert


def test_basic_poly_degree(basic_poly):
    assert basic_poly._polynom_degree([0, 1, 1, 0, 0]) == 2


@pytest.mark.parametrize('poly,factored_poly', [
    ([0, 0, 10, 5, 0, 0], [10, 5]),
    ([0, 0, 0, 0], []),
])
def test_basic_factor_out_poly(basic_poly, poly, factored_poly):
    assert basic_poly._factor_out_polynom(poly) == factored_poly


@pytest.mark.parametrize('a_range,b_range,expected_results', [
    ([[[-1, 2], [-2, 1], [1, 3]]], [[[1, 4]]], {
        (((-1, -2, 1),), ((1,),)),
        (((-1, -2, 1),), ((2,),)),
        (((-1, -2, 1),), ((3,),)),
        (((-1, -2, 2),), ((1,),)),
        (((-1, -2, 2),), ((2,),)),
        (((-1, -2, 2),), ((3,),)),
        (((-1, -1, 1),), ((1,),)),
        (((-1, -1, 1),), ((2,),)),
        (((-1, -1, 1),), ((3,),)),
        (((-1, -1, 2),), ((1,),)),
        (((-1, -1, 2),), ((2,),)),
        (((-1, -1, 2),), ((3,),)),
        (((-1, 0, 1),), ((1,),)),
        (((-1, 0, 1),), ((2,),)),
        (((-1, 0, 1),), ((3,),)),
        (((-1, 0, 2),), ((1,),)),
        (((-1, 0, 2),), ((2,),)),
        (((-1, 0, 2),), ((3,),)),
        (((0, -2, 1),), ((1,),)),
        (((0, -2, 1),), ((2,),)),
        (((0, -2, 1),), ((3,),)),
        (((0, -2, 2),), ((1,),)),
        (((0, -2, 2),), ((2,),)),
        (((0, -2, 2),), ((3,),)),
        (((0, -1, 1),), ((1,),)),
        (((0, -1, 1),), ((2,),)),
        (((0, -1, 1),), ((3,),)),
        (((0, -1, 2),), ((1,),)),
        (((0, -1, 2),), ((2,),)),
        (((0, -1, 2),), ((3,),)),
        (((0, 0, 1),), ((1,),)),
        (((0, 0, 1),), ((2,),)),
        (((0, 0, 1),), ((3,),)),
        (((0, 0, 2),), ((1,),)),
        (((0, 0, 2),), ((2,),)),
        (((0, 0, 2),), ((3,),)),
        (((1, -2, 1),), ((1,),)),
        (((1, -2, 1),), ((2,),)),
        (((1, -2, 1),), ((3,),)),
        (((1, -2, 2),), ((1,),)),
        (((1, -2, 2),), ((2,),)),
        (((1, -2, 2),), ((3,),)),
        (((1, -1, 1),), ((1,),)),
        (((1, -1, 1),), ((2,),)),
        (((1, -1, 1),), ((3,),)),
        (((1, -1, 2),), ((1,),)),
        (((1, -1, 2),), ((2,),)),
        (((1, -1, 2),), ((3,),)),
        (((1, 0, 1),), ((1,),)),
        (((1, 0, 1),), ((2,),)),
        (((1, 0, 1),), ((3,),)),
        (((1, 0, 2),), ((1,),)),
        (((1, 0, 2),), ((2,),)),
        (((1, 0, 2),), ((3,),)),
    }),
    ([[[-1, 1], [-2, -1]], [[6, 8]]], [[[1, 2]], [[4, 6]]], {
        (((-1, -2), (6,)), ((1,), (4,))),
        (((-1, -2), (6,)), ((1,), (5,))),
        (((-1, -2), (7,)), ((1,), (4,))),
        (((-1, -2), (7,)), ((1,), (5,))),
        (((0, -2), (6,)), ((1,), (4,))),
        (((0, -2), (6,)), ((1,), (5,))),
        (((0, -2), (7,)), ((1,), (4,))),
        (((0, -2), (7,)), ((1,), (5,))),
    }),
    ([[[1, 2]]], [[[4, 5], [2, 4]]], {
        (((1,),), ((4, 2),)),
        (((1,),), ((4, 3),))
    })
])
def test_basic_polys_generator(basic_poly, a_range, b_range, expected_results):
    generated_results = set([ (pa, pb) for contfrac, pa, pb in basic_poly.polys_generator(a_range, b_range)[0] ])
    assert generated_results.symmetric_difference(set(expected_results)) == set()
    is_generated_contfracs = [ isinstance(contfrac, enum_poly_params.cont_fracs.ContFrac)
                               for contfrac, pa, pb in basic_poly.polys_generator(a_range, b_range)[0] ]
    assert all(is_generated_contfracs)


@pytest.mark.parametrize('a_range,b_range,expected_results', [
    ([[[1, 3], [9, 10]]], [[[4, 5], [2, 4], [1, 2]]], {
        (((9, 10),), ((6, 8, 7),)),
        (((9, 10),), ((7, 10, 8),)),
        (((9, 11),), ((6, 8, 7),)),
        (((9, 11),), ((7, 10, 8),)),
    })
])
def test_indexed_polys_generator(a_range, b_range, expected_results):
    indexed_poly = enum_poly_params.IndexedParameterEnumPolyParams(num_of_iterations=10, enum_only_exp_conv=False,
                                                                   avoid_int_roots=True, should_gen_contfrac=True,
                                                                   avoid_zero_b=True, threshold=None, prec=80)
    generated_results = set([ (pa, pb) for contfrac, pa, pb in indexed_poly.polys_generator(a_range, b_range)[0] ])
    assert generated_results.symmetric_difference(set(expected_results)) == set()
    is_generated_contfracs = [ isinstance(contfrac, enum_poly_params.cont_fracs.ContFrac)
                               for contfrac, pa, pb in indexed_poly.polys_generator(a_range, b_range)[0] ]
    assert all(is_generated_contfracs)


@pytest.mark.parametrize('a_range,b_range,n,expected_results', [
    ([[[1, 3], [9, 10]]], [[[4, 5], [2, 4], [1, 2]]], (3, 3), {
        (((1, 9, 0),), ((4, 2, 1),)),
        (((1, 9, 0),), ((4, 3, 1),)),
        (((2, 9, 0),), ((4, 2, 1),)),
        (((2, 9, 0),), ((4, 3, 1),)),
        (((1, 0, 9),), ((4, 2, 1),)),
        (((1, 0, 9),), ((4, 3, 1),)),
        (((2, 0, 9),), ((4, 2, 1),)),
        (((2, 0, 9),), ((4, 3, 1),)),
        (((0, 1, 9),), ((4, 2, 1),)),
        (((0, 1, 9),), ((4, 3, 1),)),
        (((0, 2, 9),), ((4, 2, 1),)),
        (((0, 2, 9),), ((4, 3, 1),)),
    }),
    ([[[1, 3]], [[6, 7]]], [[[10, 11]]], (2, 3), {
        (((1, 0), (6, 0)), ((10, 0, 0),)),
        (((2, 0), (6, 0)), ((10, 0, 0),)),
        (((1, 0), (6, 0)), ((0, 10, 0),)),
        (((2, 0), (6, 0)), ((0, 10, 0),)),
        (((1, 0), (6, 0)), ((0, 0, 10),)),
        (((2, 0), (6, 0)), ((0, 0, 10),)),
        (((1, 0), (0, 6)), ((10, 0, 0),)),
        (((2, 0), (0, 6)), ((10, 0, 0),)),
        (((1, 0), (0, 6)), ((0, 10, 0),)),
        (((2, 0), (0, 6)), ((0, 10, 0),)),
        (((1, 0), (0, 6)), ((0, 0, 10),)),
        (((2, 0), (0, 6)), ((0, 0, 10),)),
        (((0, 1), (6, 0)), ((10, 0, 0),)),
        (((0, 2), (6, 0)), ((10, 0, 0),)),
        (((0, 1), (6, 0)), ((0, 10, 0),)),
        (((0, 2), (6, 0)), ((0, 10, 0),)),
        (((0, 1), (6, 0)), ((0, 0, 10),)),
        (((0, 2), (6, 0)), ((0, 0, 10),)),
        (((0, 1), (0, 6)), ((10, 0, 0),)),
        (((0, 2), (0, 6)), ((10, 0, 0),)),
        (((0, 1), (0, 6)), ((0, 10, 0),)),
        (((0, 2), (0, 6)), ((0, 10, 0),)),
        (((0, 1), (0, 6)), ((0, 0, 10),)),
        (((0, 2), (0, 6)), ((0, 0, 10),)),

    })
])
def test_sparse_polys_generator(a_range, b_range, n, expected_results):
    indexed_poly = enum_poly_params.SparseParameterEnumPolyParams(num_of_iterations=10, enum_only_exp_conv=False,
                                                                  avoid_int_roots=True, should_gen_contfrac=True,
                                                                  avoid_zero_b=True, threshold=None, prec=80,
                                                                  special_params=n)
    generated_results = set([ (pa, pb) for contfrac, pa, pb in indexed_poly.polys_generator(a_range, b_range)[0] ])
    assert generated_results.symmetric_difference(set(expected_results)) == set()
    is_generated_contfracs = [ isinstance(contfrac, enum_poly_params.cont_fracs.ContFrac)
                               for contfrac, pa, pb in indexed_poly.polys_generator(a_range, b_range)[0] ]
    assert all(is_generated_contfracs)
