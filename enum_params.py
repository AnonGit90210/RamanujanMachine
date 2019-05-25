"""Integrates most part and provides the main module: the one that implements meet-in-the-middle by exhausting the
continued fractions and creating a result hashtable, exhausting the LHS (functions applied on the target constant) and
finding 'clicks', filter unique ones and ones that converge fast enough and then refine the clicks to those that match
with the continued fractions up to an arbitrary decimal precision."""

import cont_fracs
from decimal_hashtable import DecimalHashTable
## # from decimal import Decimal as dec
from mpmath import mpf as dec, isnormal
## # from decimal import DivisionByZero, DivisionUndefined, DivisionImpossible, InvalidOperation
from gen_consts import gen_pi_const
import csv
from latex import latex_cont_frac
from postprocfuncs import EVALUATED_POSTPROC_FUNCS ,POSTPROC_FUNCS, INVERSE_POSTPROC_PAIRS ,POSTPROC_FUNCS_LATEX
from enum_poly_params import BasicEnumPolyParams, set_precision
import io
import shutil
import numpy
# The usage of progressbar adds about 1:05 minutes for a 25 bits enumeration
import progressbar
import itertools
from utils import MathOperations, grouper

# from lhs_evaluators import ULCDEnumerator, ULCDEvaluator, RationalFuncEnumerator, RationalFuncEvaluator
# converts #
# ENUMEATOR_TYPES = {'ulcd': ULCDEnumerator,
#                     'rationalfunc': RationalFuncEnumerator}
# EVALUATOR_TYPES = {'ulcd': ULCDEvaluator,
#                    'rationalfunc': RationalFuncEvaluator}
#
# EVALUATOR2TYPE = { v: k for k, v in EVALUATOR_TYPES.items() }
# ENUMERRATOR2TYPE = { v: k for k, v in ENUMERATOR_TYPES.items() }

# do we still want to keep the default values here? (some of them require special imports that are not needed otherwise)
class MITM:
    # TODO: right documentation. What does this class do?
    def __init__(self, target_generator=gen_pi_const, target_name='pi', postproc_funcs_filter=(0, 1, 2, 10),
                 trunc_integer=True, hashtable_prec=15, ab_poly_class=BasicEnumPolyParams, ab_poly_special_params=None,
                 enum_only_exp_conv=True, num_of_iterations=300, prec=50, postproc_funcs=EVALUATED_POSTPROC_FUNCS,
                 postproc_funcs_text=POSTPROC_FUNCS, postproc_funcs_text_inverse=INVERSE_POSTPROC_PAIRS,
                 postproc_funcs_latex=POSTPROC_FUNCS_LATEX):
        """target_generator - function that generates the target value.
        target_name - the name of the target, used for exporting and displaying the results.
        postproc_funcs - a list of functions that can be applied on the continued fraction results before saving to
                         the hashtable.
        postproc_funcs_filter - indices of functions from postproc_funcs to actually be applied on the resulting
                                continued fractions. The reason for this structure is that it is easier to provide
                                 backward compatibility and flexibility, as all the possible needed data (the postproc
                                 functions in this case) is stored in the class.
        trunc_integer - whether to ignore an integer difference when comparing to the target value. E.g. if
                        contfrac_val=5.2348 and target_val=3.2348, it will be considered a match (or a "click").
                        This is done in two different places - when generating the continued fractions, they're saved
                        to the hashtable without an integer-value of 0 (0.xyzaba...). When enumerating the LHS and
                        looking for clicks, they are truncated to fractional value too before being sought in the
                        hashtable.
        hashtable_prec - what decimal precision to be used in the hashtable (this includes ALL decimals, before and
                         after the decimal point).
        ab_poly_class - a class to be used for initializing the enumerator of the  a,b polynomials (of the continued
                        fraction).
        enum_only_exp_conv - whether to skip enumeration over non-exponential converging continued fractions.
        num_of_iterations - number of iterations to produce when building the contfracs for generating the hashtable.
                            This number is closely related to the hastable_prec value, as only contfracs that achieve
                            'hashtable_prec' precision in 'num_of_iterations' iterations will be stored properly and
                            may be later found.
                            In future, if we could easily estimate the rate of convergence, this parameter should be
                            replaced with an upper bound for iterations, and the required precision should be used to
                            estimate the number of required iterations.
        prec - the decimal precision to be used during calculations. This of course bounds all other used precisions."""
        self.rhs_polys_enumer = ab_poly_class(num_of_iterations=num_of_iterations,
                                              enum_only_exp_conv=enum_only_exp_conv, avoid_int_roots=True,
                                              should_gen_contfrac=False, avoid_zero_b=True, prec=prec,
                                              special_params=ab_poly_special_params)
        self._polys_enumer_num_of_iterations = num_of_iterations
        self._polys_enumer_enum_only_exp_conv = enum_only_exp_conv

        set_precision(prec)
        self.target_generator = target_generator
        self.target_name = target_name
        self.postproc_funcs_filter = postproc_funcs_filter
        self.trunc_integer = trunc_integer
        self.hashtable_prec = hashtable_prec
        self.prec = prec
        self.dec_hashtable = DecimalHashTable(self.hashtable_prec)
        self.postproc_funcs = postproc_funcs
        self.postproc_funcs_text = postproc_funcs_text
        self.postproc_funcs_text_inverse = postproc_funcs_text_inverse
        self.postproc_funcs_latex = postproc_funcs_latex
        self.filtered_params = []

    def redefine_settings(self, target_generator=gen_pi_const, target_name='pi', postproc_funcs_filter=(0, 1, 2, 10),
                          trunc_integer=True, hashtable_prec=15, ab_poly_class=BasicEnumPolyParams,
                          ab_poly_special_params=None, prec=50, postproc_funcs=EVALUATED_POSTPROC_FUNCS,
                          postproc_funcs_text=POSTPROC_FUNCS, postproc_funcs_text_inverse=INVERSE_POSTPROC_PAIRS,
                          postproc_funcs_latex=POSTPROC_FUNCS_LATEX):
        """Redefines the class's settings. See the help of '__init__' for more info.
        Two important differences are:
            1) The contfrac's a,b enumerator is recreated if ab_poly_class is changed to a
               different one.
            2) The hashtable precision is UPDATED, but the old hashtable IS NOT DELETED."""
        # TODO: document the differences between __init__ and here (ab_poly_class? hashtable updated instead of recreated).
        if not isinstance(self.rhs_polys_enumer, ab_poly_class):
            self.rhs_polys_enumer = ab_poly_class(prec=self.prec,
                                                  enum_only_exp_conv=self._polys_enumer_enum_only_exp_conv,
                                                  avoid_int_roots=True,
                                                  should_gen_contfrac=False,
                                                  num_of_iterations=self._polys_enumer_num_of_iterations,
                                                  threshold=None, special_params=ab_poly_special_params)

        set_precision(prec)
        self.target_generator = target_generator
        self.target_name = target_name
        self.postproc_funcs_filter = postproc_funcs_filter
        self.trunc_integer = trunc_integer
        self.hashtable_prec = hashtable_prec
        self.dec_hashtable.update_accuracy(self.hashtable_prec)
        self.postproc_funcs = postproc_funcs
        self.postproc_funcs_text = postproc_funcs_text
        self.postproc_funcs_text_inverse = postproc_funcs_text_inverse
        self.postproc_funcs_latex = postproc_funcs_latex
        self.filtered_params = []


    def build_hashtable(self, range_a, range_b, validate_entry_doesnt_exist=True):
        """Build a hashtable using the ranges for a, b. See the help for BasicEnumPolyParams.polys_generator for how
         range_a, range_b should be supplied.
         For each result continued fraction all the selected postproc functions will be applied on, and these results
         will be saved to the hashtable."""
        # pg - polynomial generator. The pg_len is only an (over)estimated one.
        import time
        t = time.time()
        pg, pg_len = self.rhs_polys_enumer.polys_generator(range_a=range_a, range_b=range_b)
        print('Took %f seconds to create pg, pg_len' % (time.time()-t))
        self._iter2hashtalbe(itr=pg, itr_len=pg_len, validate_entry_doesnt_exist=validate_entry_doesnt_exist)

    def _iter2hashtalbe(self, itr, itr_len, validate_entry_doesnt_exist=True, print_problematic=False):
        """Converts the (a, b) polynomials generator itr of an (estimated) length itr_len to a contfrac hashtable."""
        filtered_postproc_funcs = { (i, self.postproc_funcs[i]) for i in self.postproc_funcs_filter }

        ab_blacklist = {}
        if validate_entry_doesnt_exist:
            for k, v in self.dec_hashtable.items():
                if k == 'parameters':
                    continue
                for ab, postproc_func_ind in v:
                    ab_blacklist.setdefault(ab, set()).add(postproc_func_ind)


        print('The following is a rough estimation only, and may be wrong to an order of 2-5 times.')
        # This was an attempt to update the progressbar only once in 50 iterations, to speed things up.
        # It turned out to slow things down :(
        # for iter_group in progressbar.progressbar(grouper(itr, 50), max_value=int(itr_len/50), poll_interval=1):
        #     for iter_element in iter_group:
        #         if iter_element is None:
        #             continue
        #         cont_frac, pa, pb = iter_element
        for pa, pb in progressbar.progressbar(itr, max_value=itr_len):
            if (pa, pb) in ab_blacklist:
                postproc_funcs_to_eval = filtered_postproc_funcs - ab_blacklist[(pa, pb)]
                # if no new postproc functions to evaluate for this pa, pb, skip
                if not postproc_funcs_to_eval:
                    continue
            else:
                postproc_funcs_to_eval = filtered_postproc_funcs

            try:
                cur_cont_frac_val = cont_fracs.eval_dec_contfrac_by_polys(pa, pb, self._polys_enumer_num_of_iterations)
            ## # except (DivisionByZero, DivisionUndefined, DivisionImpossible, InvalidOperation):
            except ZeroDivisionError:
                continue
            ## # # is_normal = finite nonzero, not NaN, with exponent > Emin=-999999999999999999
            ## # if not cur_cont_frac_val.is_normal():
            # is_normal = finite nonzero, not NaN
            if not isnormal(cur_cont_frac_val):
                if print_problematic:
                    print('problematic number')
                    print(cur_cont_frac_val)
                continue

            for post_func_ind, post_f in postproc_funcs_to_eval:
                # print(cur_cont_frac_val)
                # Seems redundant - we're already enumerating only filtered functions
                # if post_func_ind not in self.postproc_funcs_filter:
                #     pass
                try:
                    k = post_f(cur_cont_frac_val)
                except:
                    print(k)
                    raise
                ## # if not k.is_normal():
                if not isnormal(k):
                    continue
                k = abs(k)
                if self.trunc_integer:
                    k -= int(k)
                if ((pa, pb), post_func_ind) not in self.dec_hashtable.get(k, []):
                    self.dec_hashtable.setdefault(k, []).append(((pa, pb), post_func_ind))
        print()

    def find_clicks(self, lhs_classes, lhs_enumerator_params):
        """Enumerating the LHS and looking for clicks. The comparison is made up to the accuracy of the stored values in
        the hashtable.
        Parameters:
        lhs_classes - a dictionary with 'enum' and 'eval' keys, pointing to the enumeration and evaluation classes in
                      lhs_evaluators.py.
        lhs_enumerator_params - parameters to feed the lhs_enumerator. See the help for the relevant enumerator for more
                                info."""
        filtered_params = []
        # lhs_enumerator = ENUMERATOR_TYPES[lhs_type](lhs_enumerator_params, self.target_generator())
        # lhs_evaluator = EVALUATOR_TYPES[lhs_type]
        lhs_enumerator = lhs_classes['enum'](lhs_enumerator_params, self.target_generator())
        lhs_evaluator = lhs_classes['eval']
        lhs_generator = lhs_enumerator.generator()
        for enum_res_obj in progressbar.progressbar(lhs_generator, max_value=len(lhs_enumerator)):
            r = abs(enum_res_obj.get_val())
            if self.trunc_integer:
                r -= int(r)
            if r in self.dec_hashtable:
                # if the lhs is degenerate (i.e. (3x+12x^4) / (x+4x^4)) then skip.
                # this is checked here, because this check may require a major additional runtime, and checking it only
                # for passing clicks is saving time (both in compare to validate in advance or not validating at all).
                if enum_res_obj.is_degenerate():
                    continue
                if r in self.dec_hashtable:
                    signed_r = r
                else:
                    signed_r = -r
                # the filtered_params' elements structure is:
                # ((a_poly, b_poly), lhs_eval_instance, post_func_ind, convergence_info)
                filtered_params.extend([ (ab, lhs_evaluator(enum_res_obj), post_func_ind, None)
                                         for ab, post_func_ind in self.dec_hashtable[signed_r] ])
        self.filtered_params = filtered_params

    def refine_clicks_with_const_num_of_iters(self, accuracy=10, num_of_iterations=3000, print_clicks=False):
        """Filters the clicks and refines them by validation up to an arbitrary accuracy.
        Parameters:
            accuracy - how many digital digits to compare (digits only. +1 is added automatically for the decimal dots).
            num_of_iterations - number of iterations (depth of the contfrac to calculate)
            print_click - print every new click that passes. If you know what's good, you'll probably keep it False."""
        contfrac_evaluator = lambda ab, accuracy, convergence_info: \
            cont_fracs.eval_dec_contfrac_by_polys(ab[0], ab[1], num_of_iterations)
        self._refine_clicks(accuracy=accuracy, contfrac_evaluator=contfrac_evaluator, print_clicks=print_clicks)

    def refine_clicks_with_convergence_info(self, accuracy=10, num_of_iterations=3000, print_clicks=False):
        """Filters the clicks and refines them by validation up to an arbitrary accuracy.
        Should be called only after filter_clicks_by_approach_type.
        Parameters:
            accuracy - how many digital digits to compare (digits only. +1 is added automatically for the decimal dots).
            num_of_iterations - upper bound for the number of iterations
            print_click - print every new click that passes. If you know what's good, you'll probably keep it False."""
        def evaluate_contfrac(ab, accuracy, convergence_info):
            cont_frac = cont_fracs.ContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            if convergence_info:
                cont_frac.set_approach_type_and_params(convergence_info)
            # cont_frac.reinitialize(a_coeffs=ab[0], b_coeffs=ab[1])
            # This is a horrible hack, slowing down EVERYTHING. As a first step, this should be replaced by a small
            # correction to the convergence parameters. Maybe taking the upper/lower bound of confidence from the
            # fitting results.
            # TODO: fix accoridng to the above comment.
            cont_frac.gen_iterations(num_of_iterations, dec('1E-%d' % (accuracy+100)))
            return cont_frac.get_result()
        self._refine_clicks(accuracy=accuracy, contfrac_evaluator=evaluate_contfrac, print_clicks=print_clicks)

    def _refine_clicks(self, accuracy, contfrac_evaluator, print_clicks):
        """Filters the clicks and refines them by validation up to an arbitrary accuracy.
        Parameters:
            accuracy - how many digital digits to compare (digits only. +1 is added automatically for the decimal dots).
            contfrac_evaluator - a function the takes ab(=[a, b]), accuracy, convergence_info and return the appropriate
                                 value of the continued fraction.
            print_click - print every new click that passes. If you know what's good, you'll probably keep it False."""
        refined_params = []
        target_value = self.target_generator()
        # cont_frac = cont_fracs.ContFrac()
        for ab, lhs_res_obj, post_func_ind, convergence_info in progressbar.progressbar(self.filtered_params):
            try:
                signed_rhs = self.postproc_funcs[post_func_ind](contfrac_evaluator(ab, accuracy, convergence_info))
            ## # except DivisionByZero:
            except ZeroDivisionError:
                ## # print('DivisionByZero exception', ab)
                print('ZeroDivisionError exception', ab)
                continue
            rhs = abs(signed_rhs)
            ## # if not rhs.is_normal():
            if not isnormal(rhs):
                continue
            signed_lhs = lhs_res_obj.get_val()
            lhs = abs(signed_lhs)
            if self.trunc_integer:
                rhs -= int(rhs)
                lhs -= int(lhs)
            if self.compare_dec_with_accuracy(rhs, lhs, accuracy):
                if print_clicks:
                    print(ab)
                    print(lhs_res_obj)
                    print(rhs)
                    print(lhs)
                    print('')
                if (signed_lhs > 0 and signed_rhs < 0) or (signed_lhs < 0 and signed_rhs > 0):
                    # print('flipping sign')
                    lhs_res_obj.flip_sign()
                    signed_lhs *= -1
                try:
                    lhs_res_obj.add_int(int(signed_rhs - signed_lhs))
                except:
                    raise RuntimeError('signed_rhs=%s, signed_lhs=%s, int(signed_rhs-signed_lhs)=%s' %
                                       (signed_rhs, signed_lhs, int(signed_rhs-signed_lhs)))
                # Turns the lhs to a canonic form. In example for a rational func, this will be the reduced form:
                # p' = p/gcd(p, q)  ,  q' = q/gcd(p, q)
                lhs_res_obj.canonalize_params()
                refined_params.append((ab, lhs_res_obj, post_func_ind, convergence_info))
        self.filtered_params = refined_params

    def filter_uniq_params(self):
        """Filters out non-unique pairs of match lhs,rhs.
        In example: if lhs_1=p(x)/q(x), rhs_1=contfrac(a(i), b(i)) and lhs_2=-p'(x)/q'(x), rhs_2=contfrac(-a(i), b(i))
        where p'/q' = p/q = (p/gcd(p,q)) / (q/gcd(p,q)).
        """
        unique_params_filtered = []
        # params = ab, lhs_res_obj, post_func_ind, convergence_info
        for params in progressbar.progressbar(self.filtered_params):
            is_unique = True
            for uniqe_params in unique_params_filtered:
                uniqe_lhs_res_obj= uniqe_params[1]
                if uniqe_lhs_res_obj.is_equiv(uniqe_params, params):
                    is_unique = False
                    break
            if is_unique:
                unique_params_filtered.append(params)
        self.filtered_params = unique_params_filtered

    def filter_pos_integer_roots_numerators(self):
        """Filters out contfracs for which any of the b interlaces have an integer positive root. The roots are found
        numerically and are considered real if Im(root) < 0.001 and integer if fractional_part(root) < 0.001"""
        valid_params = []
        for params in self.filtered_params:
            b = params[0][1]
            roots = itertools.chain([ r for b_p in b for r in numpy.roots(b_p[::-1]) ])
            real_int_roots = [ ((abs(r.imag) < 0.001) and (r.real - r.real.round() < 0.001) and (r.real > 0.001))
                               for r in roots ]
            if not any(real_int_roots):
                valid_params.append(params)
        self.filtered_params = valid_params

    def filter_only_exp_convergence(self, print_surprising_nonexp_contfracs=False):
        """Filters out only contfracs that converge (at least) exponentially.
        Parameters:
            print_surprising_nonexp_contfracs - True for printing contfracs that DO NO converge exponentially."""
        # TODO: add either here or in filter_clicks_by_approach_type the discriminant test for exp. convergence
        params_list = self.filtered_params
        filtered_params_list = []
        for cf_params in progressbar.progressbar(params_list):
            ab, lhs_res_obj, post_func_ind, convergence_info = cf_params
            cont_frac = cont_fracs.ContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            if cont_frac.is_convergence_fast():
                filtered_params_list.append(cf_params)
            elif print_surprising_nonexp_contfracs:
                print('Surprising non-exponential convergence continuous fraction:')
                print(cf_params)
                cont_frac.estimate_approach_type_and_params()
                print(cont_frac.get_approach_type_and_params())
        self.filtered_params = filtered_params_list

    def filter_clicks_by_approach_type(self, whitelist=['exp', 'super_exp', 'fast'], blacklist=None):     # , filter_uniq_list=True):
        """Filters only clicks withe convergence from the whitelist, or drops clicks with convergence from the black
         list. One of whitelist/black is required, and one only. See help for
         ContFrac._estimate_approach_type_and_params_inner_alg
         for more info about convergence types."""
        if not any([whitelist, blacklist]):
            raise ValueError('One is required: whitelist, blacklist')
        if whitelist and blacklist:
            raise ValueError('Only one is possible: whitelist, blacklist')
        params_list = self.filtered_params

        filtered_params_list = []
        for cf_params in progressbar.progressbar(params_list):
            ab, lhs_res_obj, post_func_ind, convergence_info = cf_params
            cont_frac = cont_fracs.ContFrac(a_coeffs=ab[0], b_coeffs=ab[1])
            try:
                cont_frac.estimate_approach_type_and_params()
                approach_type, approach_params = cont_frac.get_approach_type_and_params()
            except Exception as e:
                print('Problems while estimating the following cf_params in "filter_clicks_by_approach_type"')
                print('Exception has occurred. Skipping.')
                print(cf_params)
                print(e)
                continue
            if (whitelist and approach_type in whitelist) or (blacklist and approach_type not in blacklist):
                cf_params = (ab, lhs_res_obj, post_func_ind, (approach_type, approach_params))
                filtered_params_list.append(cf_params)
        self.filtered_params = filtered_params_list

    def get_filtered_params(self):
        """Returns the filtered parameters (the found 'clicks')."""
        return self.filtered_params

    def reset_hashtable(self):
        """Deletes the hashtable and creates a new one (with the earlier, predefined precision)."""
        del self.dec_hashtable
        self.dec_hashtable = DecimalHashTable(self.hashtable_prec)

    def get_results_as_eqns(self, postproc_funcs, ignore_zerob_exceptions=True, depth=5):
        """Exports the results to a list of latex equations.
        postproc_funcs - a list of (textual) postproc funcs. Needs to be given as a parameter since so far, only the
                         evaluated postproc funcs were supplied. This may be replaced by adding to postprocfuncs.py
                         beside POSTPROC_FUNCS_LATEX a similar dictionary from the evaluated functions to the latex
                         generators.
        ignore_zerob_exceptions - if a ZeroB exception is raised while building the contfracs, ignore it and continue
                                  to the next contfrac. This actually shouldn't happen, as these contfracs shouldn't be
                                  saved in the first place.
        depth - to what depth should the contfrac expression be built.
        Returns: """
        eqns = []
        eval_poly = MathOperations.subs_in_polynom
        # target names that should be replaces by anothers when exporting
        known_targets = {'pi': '\pi',
                        'phi': r'\varphi'}

        for ab, lhs_res_obj, post_func_ind, convergence_info in self.filtered_params:
            pa, pb = ab
            try:
                cont_frac, postproc_res, lhs_res_obj = self.build_contfrac_from_params((ab, lhs_res_obj, post_func_ind,
                                                                                        convergence_info))
            except cont_fracs.ZeroB:
                if ignore_zerob_exceptions:
                    continue
                else:
                    raise

            # evaluate the first 'depth' values of a and b
            a = [eval_poly(pa[i % len(pa)], i) for i in range(depth)]
            b = [eval_poly(pb[i % len(pb)], i) for i in range(depth)]

            # rename the target name, if needed (e.g. 'pi'-->'\pi')
            if self.target_name in known_targets:
                target_name = known_targets[self.target_name]
            else:
                target_name = self.target_name

            # Creates the equation object. At first lhs,rhs are latex expressions, then they are changed according to
            # the postproc function.
            lhs = lhs_res_obj.get_latex_exp(target_name)
            rhs = latex_cont_frac(a, b)
            lhs, rhs = self.postproc_funcs_latex[postproc_funcs[post_func_ind]](lhs, rhs)
            eqn = r'{0} = {1}'.format(lhs, rhs)

            # Appends equation
            eqns.append(eqn)

        return eqns

    def export_to_csv(self, filename):
        """Exports the found results to a csv file."""
        with io.StringIO(newline='') as csvbuffer:
            csvwriter = csv.writer(csvbuffer)
            # csvwriter.writerow(['postproc_funcs', postfuncs, 'target_name', self.target_name])
            # First row: "target_name; 'pi'". The above line can be uncommented to save the postproc funcs too. However,
            # if the code is used correctly (i.e. postproc funcs are only added, never, moved around or deleted), then
            # the later specified postproc_func indices should be enough.
            csvwriter.writerow(['target_name', self.target_name])
            # Second row: columns titles for all the following results
            csvwriter.writerow(['a poly2sympoly [a_0, a_1, ...]', 'b poly2sympoly  [b_0, b_1, ...]', 'postproc_func',
                                'LHS type', 'LHS params',
                                'convergence type', 'convergence rate', 'postfunc(cont_frac)',
                                'LHS val'])
            # add a row for each result
            for ab, lhs_res_obj, post_func_ind, convergence_info in self.filtered_params:
                pa, pb = ab
                # lhs_type = EVALUATOR2TYPE[type(lhs_res_obj)]
                lhs_type =str(lhs_res_obj)
                # builds the contfrac from the results to receive the final, numerical result (after the postproc func).
                cont_frac, postproc_res, lhs_res_obj = self.build_contfrac_from_params((ab, lhs_res_obj,
                                                                                        post_func_ind,
                                                                                        convergence_info))
                # csvwriter.writerow([pa, pb, post_func_ind, lhs_type, lhs_res_obj.get_params(),
                #                     convergence_info[0], convergence_info[1],
                #                     postproc_res.to_eng_string(), lhs_res_obj.get_val().to_eng_string()])
                csvwriter.writerow([pa, pb, self.postproc_funcs[post_func_ind], lhs_type, lhs_res_obj.get_params(),
                                    convergence_info[0], convergence_info[1],
                                    str(postproc_res), str(lhs_res_obj.get_val())])
            # Saves the results
            with open(filename, 'w', newline='') as csvfile:
                csvbuffer.seek(0)
                shutil.copyfileobj(csvbuffer, csvfile)

    def build_contfrac_from_params(self, params, iterations=400):
        """Builds a continued fraction, post-processed value and an LHS object from params.
           params - (ab, lhs_res_obj, post_func_ind, convergence_info)
                    ab - [a_poly, b_poly] where a_poly/b_poly is of the format [poly1, poly2, ...] for these polynomials
                         to be used as an interlace, and poly1/2/... is of the format
                         [c0, c1, c2, ...] <=> c0+c1*x+c2*x^2+...
                    lhs_res_obj - an LHS object
                    post_func_ind - index of a postproc func in postproc_funcs list to be use.
                    convergence_info - convergence info produced by the original contfrac object. UNUSED.
           Returns: (cont_frac, postproc_func(cont_frac_value), lhs_res_obj)
          """
        ab, lhs_res_obj, post_func_ind, convergence_info = params
        pa, pb = ab
        cont_frac = cont_fracs.ContFrac(a_coeffs=pa, b_coeffs=pb)
        cont_frac.gen_iterations(iterations)
        return cont_frac, self.postproc_funcs[post_func_ind](cont_frac.get_result()), lhs_res_obj

    @staticmethod
    def compare_dec_with_accuracy(d1, d2, accuracy):
        """Compares two decimals up to a specified accuracy."""
        # Decimal.quantize('1.0000...') may be used too. Profiling shows, however, that string a conversion and slicing
        # is faster.
        # +1 for decimal dot
        accuracy += 1
        return str(d1)[:accuracy+1] == str(d2)[:accuracy+1]

