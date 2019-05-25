"""Implementation of enumeration over different polynomials types, for the a, b polynomials of a continued fractions."""

import itertools
from functools import wraps
import cont_fracs
## # from decimal import Decimal as dec
## # import decimal
from mpmath import mpf as dec, isint
from mpmath import isint
import mpmath
from scipy.special import binom

def set_precision(prec):
    """Sets the precision of the current decimal context"""
    ## # decimal.getcontext().prec=prec
    mpmath.mp.dps = prec

def _len_decorator(func):
    """Used as a decorator to estimated the length of a generator with signature:
    gen(range_a, range_b, prec)
    Instead of returning a new generator 'gen',
    the tuple (gen, estimated_len(gen)) is returned."""
    @wraps(func)
    def wrapper(self, range_a, range_b, prec=None):
        gen_len = 1
        for ar in range_a:
            for r in ar:
                gen_len *= (r[1] - r[0])
        for br in range_b:
            for r in br:
                gen_len *= (r[1] - r[0])
        if hasattr(self, 'len_factor_func'):
            return func(self, range_a, range_b, prec), gen_len * self.len_factor_func(range_a, range_b)
        else:
            return func(self, range_a, range_b, prec), gen_len
    return wrapper


class NormalMetaClass(type):
    """Provides the value of str(BasicEnumPolyParams) => 'normal'"""
    def __str__(self):
        return 'normal'


class BasicEnumPolyParams(metaclass=NormalMetaClass):
    """Enumerates over a,b polynomials for continued fractions generation. The generated polynomials are taken as is,
    unlike the non-basic enum poly2sympoly params.
    Other classes may inherit from this one and implement the static method manipulate_poly as a generator, to generate
    from each "regular polynom" the required polynomials (or yield a single transformed polynomial)."""
    def __init__(self, num_of_iterations=300, enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True,
                 avoid_zero_b=True, threshold=None, prec=80, special_params=None):
        """
        init
        :param self: self
        :param num_of_iterations: if should_gen_contfrac=True, then how many iterations should be promoted.
        :param enum_only_exp_conv: skip contfracs that surely won't converge (at least) exponentially
        :param avoid_int_roots: avoid contfracs for which b (or any of its interlaces) have an integer root.
                                it will affect only b polynomial of degree<=3.
                                notice that if interlace is used, this may skip valid contfracs (the i root value may be
                                assigned to another interlace polynomial instead of the one with the root).
        :param should_gen_contfrac: generate ContFrac object and return (cont_frac, pas, pbs) instead of (pas, pbs)
        :param avoid_zero_b: raise exception ZeroB and cancel if finite contfrac (b_i=0 for some i) is achieved
        :param threshold: threshold for considering contfrac_res == target_val
                          if abs(contfrac_res-target_val) < threshold.
                          Defaults to 10^-4.
        :param prec: decimal precision to be used for calculations
        :param special_params: unused. Used differently by inheriting subclasses.
        :return: nothing
        """
        self.num_of_iterations = num_of_iterations
        self._enum_only_exp_conv = enum_only_exp_conv
        self._avoid_int_roots = avoid_int_roots
        self._should_gen_contfrac = should_gen_contfrac
        self._avoid_zero_b = avoid_zero_b
        self.good_params = []
        if threshold is None:
            self.threshold = dec(10)**dec(-4)
        else:
            self.threshold = threshold
        self._prec = prec
        set_precision(prec)

    def reset_precision(self):
        """Resets the used decimal precision to the predefined (during initialization) precision."""
        set_precision(self._prec)

    def reinitialize_good_params(self):
        """Reinitializes the good parameters (the sets of a,b that were found and saved by the 'enum_params' method)."""
        self.good_params = []

    @_len_decorator
    def polys_generator(self, range_a, range_b, prec=None):
        """Generates the polynomials from range_a and range_b.
        Parameters:
            range_a - for example: [ [[], []], [[], [], []], [[], [], [], []], [[], []] ] is a 4-interlace
                      with degrees of 2,3,4,2 . [m n] means running on coefficients between m to n-1.
            range_b - as range_a.
            prec - ignored. Left for compatiblity with the _len_decorator decoration.
        Returns:
            yields pairs for (a_poly, b_poly), each of them is a list of polynomials (for interlace)."""
        # Create an instance. To improve runtime, this method will be reinitialized over and over instead of creating.
        # new instances every so often.
        cont_frac = cont_fracs.ContFrac([1], [1], avoid_zero_b=self._avoid_zero_b)
        # TODO: delete the following line. It left as a documentation backup until everything works.
        # cont_frac = cont_fracs.ContFrac([1]*range_a, [0]*range_b, avoid_zero_b=self._avoid_zero_b)

        # Catersian product of all the possibilities of a.
        a_params_iterator = itertools.product(*[ itertools.product(*[ range(*r) for r in ra ]) for ra in range_a ])
        for pas_premanipulate in a_params_iterator:
            # In more complex cases, manipulated versions of the a polynomial may be wished.
            # Create a generator for this manipulated versions.
            pas_manipulated_gen = self.manipulate_poly(pas_premanipulate, 'a')
            for pas in pas_manipulated_gen:
                # Cartesian product of b.
                b_params_iterator = itertools.product(*[ itertools.product(*[ range(*r) for r in rb ])
                                                         for rb in range_b ])
                for pbs_premanipulate in b_params_iterator:
                    # Manipulated versions for b.
                    pbs_manipulated_gen = self.manipulate_poly(pbs_premanipulate, 'b')
                    for pbs in pbs_manipulated_gen:
                        # Check for integers roots of b polynomials and avoid if required.
                        if self._avoid_int_roots and any([self._does_have_int_positive_roots(pb) for pb in pbs]):
                            continue
                        # in the case of an interlace in which all the interlace-polynomials are identical,
                        # squeeze it to a single polynomial with no interlace
                        if len(pas) > 1 and all([ pas[0] == p for p in pas ]):
                            pas = (pas[0],)
                        if len(pbs) > 1 and all([ pbs[0] == p for p in pbs ]):
                            pbs = (pbs[0],)
                        # if no interlace and only exponential convergence contfracs should be enumerated, make sure
                        # that 2*deg(a) >= deg(b)
                        if len(pas) == 1 and len(pbs) == 1 and self._enum_only_exp_conv and \
                           self._polynom_degree(pbs[0]) > 2 * self._polynom_degree(pas[0]):
                            continue
                        # generate contfrac and return everything / return the polynomials
                        if self._should_gen_contfrac:
                            cont_frac.reinitialize(pas, pbs)
                            try:
                                cont_frac.gen_iterations(self.num_of_iterations)
                            except cont_fracs.ZeroB:
                                continue
                            yield (cont_frac, pas, pbs)
                        else:
                            yield (pas, pbs)

    @staticmethod
    def _does_have_int_positive_roots(poly):
        """For a poly2sympoly of deg(poly2sympoly) < 3 (after factoring out any power of x^n), check if it has integer
        positive roots.
        Returns a boolean."""
        poly = BasicEnumPolyParams._factor_out_polynom(poly)
        if len(poly) == 0:
            return True
        elif len(poly) == 1:
            return False
        ## # elif len(poly) == 2 and (poly[0] / poly[1]).is_integer() and (poly[0] > 0) != (poly[1] > 0):
        elif len(poly) == 2 and isint(poly[0] / poly[1]) and (poly[0] > 0) != (poly[1] > 0):
            return True
        elif len(poly) == 3:
            discrim = (poly[1]**2-4*poly[0]*poly[2])
            if discrim < 0:
                return False
            discrim **= 0.5
            ## # if ((discrim - poly[1])/(2*poly[0])).is_integer():
            if isint((discrim - poly[1])/(2*poly[0])) and ((discrim - poly[1])/(2*poly[0]) > 0):
                return True
            if isint((-discrim - poly[1])/(2*poly[0])) and ((-discrim - poly[1])/(2*poly[0]) > 0):
                return True

        return False

    @staticmethod
    def _polynom_degree(p):
        """Finds deg(p). It may be different than len(p). E.g. deg([1, 1, 0, 0]) is 1 and not 3."""
        deg = len(p)
        for i in p[::-1]:
            if i != 0:
                break
            deg -= 1
        return deg - 1

    @staticmethod
    def _factor_out_polynom(poly):
        """In: poly2sympoly - a polynomial as a list of coefficients
        Out: a polynomial factored out of x^n. In example: [0, 0, 1, 2, 0, 5, 0] --> [1, 2, 0, 5]"""
        return list(itertools.takewhile(lambda x: x != 0, itertools.dropwhile(lambda x: x == 0, poly)))

    @staticmethod
    def manipulate_poly(poly, is_a_or_b_poly):
        yield poly


class IndexedMetaClass(NormalMetaClass):
    """Provides the value of str(IndexedParameterEnumPolyParams) => 'indexed'"""
    def __str__(self):
        return 'indexed'


class IndexedParameterEnumPolyParams(BasicEnumPolyParams, metaclass=IndexedMetaClass):
    """This class provides enumeration over polynomials of the pattern:
        a_0*x^n + a_1*(x+1)^n + a_2*(x+2)^n + a_n*(x+n)^n
     The input polynomial parameters are considered as a_0,...,a_n of the above pattern."""
    def __init__(self, num_of_iterations=300, enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True,
                 avoid_zero_b=True, threshold=None, prec=80, special_params=None):
        super().__init__(num_of_iterations=num_of_iterations, enum_only_exp_conv=enum_only_exp_conv,
                         avoid_int_roots=avoid_int_roots, should_gen_contfrac=should_gen_contfrac,
                         avoid_zero_b=avoid_zero_b, threshold=threshold, prec=prec, special_params=special_params)

    @staticmethod
    def manipulate_poly(poly, is_a_or_b_poly):
        poly_new = []
        # run over the different interlace polynomials
        for p in poly:
            p_new = [ 0 ] * len(p)
            deg_p = len(p)-1
            # p_new += sum_i a_i*(x+i)**deg_p
            # for each a_i, calculate the the contribution of a_i*(x+i)^n to the coefficients of the expanded polynomial
            for i, a_i in enumerate(p):
                if a_i == 0:
                    continue
                # p_new += (x+i)**deg_p
                for j in range(deg_p+1):
                    p_new[j] += a_i * int(binom(deg_p, j)) * i**(deg_p-j)
            poly_new.append(tuple(p_new))
        yield tuple(poly_new)


class SparseMetaClass(NormalMetaClass):
    """Provides the value of str(SparseParameterEnumPolyParams) => 'sparse'"""
    def __str__(self):
        return 'sparse'


class SparseParameterEnumPolyParams(BasicEnumPolyParams, metaclass=SparseMetaClass):
    """This class provides enumeration over sparse polynomials, such as: 5*x^2+4*x^15, 1+x^7+3*x^20"""
    def __init__(self, num_of_iterations=300, enum_only_exp_conv=False, avoid_int_roots=True, should_gen_contfrac=True,
                 avoid_zero_b=True, threshold=None, prec=80, special_params=None):
        """special_params - (n_a, n_b), where n is the maximum degree + 1 (the maximum "size") of the polynomial,
        and k is the number of nonzero elements in the sparse polynomial which is assumed to be the number of enumerated
        coefficients of the first interlace polynomial. All interlaces are assumed to have equal number of coefficients!
        E.g. for special_params=(2, 2_ with coefficients of [[1, 3]] for a, we get for a: 1, 2, x, 2x, x^2, 2x^2"""
        # for the calculation of the sparse polynom by "n over k" options
        self.n_a, self.n_b = special_params
        # give a product factor of the length, for given ranges of a, b
        self.len_factor_func = lambda range_a, range_b: int(binom(self.n_a, len(range_a[0])) *
                                                            binom(self.n_b, len(range_b[0])))
        super().__init__(num_of_iterations=num_of_iterations, enum_only_exp_conv=enum_only_exp_conv,
                         avoid_int_roots=avoid_int_roots, should_gen_contfrac=should_gen_contfrac,
                         avoid_zero_b=avoid_zero_b, threshold=threshold, prec=prec, special_params=special_params)

    def manipulate_poly(self, poly, is_a_or_b_poly):
        # creates a set of [j_1, j_2, ..., j_k] with j_i in range(n) for all i's. E.g. [2, 4, 5] for n=6, j=3
        k = len(poly[0])
        n = self.n_a if is_a_or_b_poly == 'a' else self.n_b
        single_poly_masks = itertools.combinations(range(n), k)
        # creates combinations of len(poly2sympoly) possible masks, one for each interlace polynomial in poly2sympoly. Sorted.
        all_polys_masks_sorted = itertools.combinations_with_replacement(single_poly_masks, len(poly))
        # permuted order of the sorted masks combinations, each permuted set of masks as a set.
        # Identical orders are compressed to a single instance.
        # i.e. for masks combination of [[1, 4], [1, 4]] this is the only permutation,
        # thus it results with {[[1, 4], [1, 4]]}
        all_polys_masks_groups = map(lambda masks: set(itertools.permutations(masks, len(poly))),
                                     all_polys_masks_sorted)

        for polys_masks_group in all_polys_masks_groups:
            for polys_masks in polys_masks_group:
                yield self._apply_poly_template_on_coeffs_list(poly, polys_masks, n)

    def _apply_poly_template_on_coeffs_list(self, coeffs_list, poly_masks, n):
        # for each interlace polynomial coefficients and its mask, create a sparse polynomial with the corresponding
        # coefficients and the right corresponding mask.
        return_value = ( tuple(( (p[mask.index(i)] if i in mask else 0) for i in range(n) ))
                 for p, mask in zip(coeffs_list, poly_masks) )
        return tuple(return_value)

# apery a: 2 + n (2 + n) (4 + 3 n) = 2 + 8 n + 10 n^2 + 3 n^3
# apery b: 2 n^5 (-1 + 2 n) = -2 n^5 + 4 n^6
# contfrac(a,b) = 5/(2apery)


def define_polys_types():
    """Autogenerates the global AB_POLYS_TYPE = {<poly_type_name>: poly_class}"""
    import sys
    from inspect import getmembers, isclass
    from itertools import groupby
    global AB_POLYS_TYPES
    polys_grouped = groupby([
                                poly_class[1] for poly_class in
                                getmembers(sys.modules[__name__],
                                           lambda cl: isclass(cl) and issubclass(cl,
                                                                                 BasicEnumPolyParams))
                                ],
                            lambda cl: str(cl))
    AB_POLYS_TYPES = {}
    for k, classes in polys_grouped:
        for c in classes:
            AB_POLYS_TYPES[k] = c

# Creates the global AB_POLYS_TYPES
AB_POLYS_TYPES = {}
define_polys_types()
