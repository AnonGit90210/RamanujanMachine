"""Different utilities with no clear belonging."""

import time
import sys
from itertools import zip_longest

class MeasureRuntime():
    """"Stopper" class for timing execution. Example:
    m = MeasureRuntime()
    m.start(measure)
    # <code_to_measure1>
    run_time1 = m.measure_time()
    # <code_to_measure2>
    run_time2 = m.measure_time()"""

    def __init__(self):
        self.current_index = 0
        self.previous_index = 0
        self.times = [0, 0]
        self.started = False

    def start_measure(self):
        """Starts the measure."""
        if self.started:
            print('Measure has already began.')
            return
        self._start_time = time.time()
        self.times[self.current_index] = self._start_time
        self.current_index = (self.current_index + 1) % 2
        self.started = True

    def measure_time(self):
        """MeasureS the time since the last call to measure_time (or start_measure)."""
        self.times[self.current_index] = time.time()
        self.previous_index = self.current_index
        self.current_index = (self.current_index + 1) % 2
        return self.times[self.previous_index] - self.times[self.current_index]

    def get_total_runtime(self):
        """Returns the total runtime since start_measure was called."""
        return time.time() - self._start_time

    def get_last_measured_time(self):
        """Returns the last measured time."""
        return self.times[self.previous_index] - self.times[self.current_index]

    def is_started(self):
        """Returns True if time is already being measured."""
        return self.started


class MathOperations:
    """Collects different math operations needed in several places, either operations that aren't supplied with Decimal
    or operations related to polynomials etc."""
    @classmethod
    def gcd(cls, a, b):
        """Calculates the Greatest Common Divisor of a and b. Useful for Decimal that has no built-in gcd.

        Unless b==0, the result will have the same sign as b (so that when
        b is divided by it, the result comes out positive).
        (Implemented here to be used with Decimal/other types that math doesn't support.)
        """
        if b > a:
            return cls.gcd(b, a)
        while b:
            a, b = b, a % b
        return a

    @staticmethod
    def subs_in_polynom(coeffs, x):
        """Substitues x in the polynomial represented by the list of coeffs
        Parameters:
            coeffs - polynoms coefficients, coeffs[0]+coeffs[1]*x+...
            x - a value to substitute"""
        return sum([coeffs[j] * x ** j for j in range(len(coeffs))])


class ProgressBar(enumerate):
    """Deprecated. Use progressbar2 instead.

    Provides an iterator with a progress bar and optionally an enumerator."""
    def __new__(cls, *args, **kwargs):
        """Creates a new instance of an iterator that prints automatically a progress bar.
        Optional parameters are:
            est_len - estimated length of the supplied generator/iterator, overwrites len(itr) if exists.
            enumerate - whether i, v should be return as by enumerate(itr), or only v. Default: False.
            update_freq - how frequently should the progress bar be updated, in seconds. May take longer if a single
                          iteration takes longer then update_freq. Default: 1."""
        should_enum = kwargs.pop('enumerate', False)
        update_frequency = kwargs.pop('update_freq', 1)
        try:
            estimated_length = kwargs.pop('est_len')
        except KeyError as exc:
            estimated_length = None

        self = super().__new__(cls, *args, **kwargs)
        if estimated_length is None:
            if not hasattr(args[0], '__len__'):
                raise ValueError("est_len=... must be supplied to ProgressBar(itr) if itr.__len__ doesn't exist.") \
                    from exc
            estimated_length = len(args[0])

        self._should_enum = should_enum
        self._update_frequency = update_frequency
        self._estimated_length = estimated_length
        # Used to follow runtime and update estimated finish time
        self._time_stopper = MeasureRuntime()
        # Used to remember how much time has passed since last estimated finish time was printed
        self._printing_time_counter = 0
        # how much time left, and how many iterations is it based on.
        self._left_time_estimation = {'left_time': 0, 'num_of_iters': 0, 'passed_time': 0}
        self._progressbar_template = '\r{0:>3}%\t|\tEstimated time left: {1:>02}:{2:>02}:{3:>02}:{4:>02}\t|\tTime passed: {5:>02}:{6:>02}:{7:>02}:{8:>02}'
        return self

    def __next__(self):
        """If enumeration is requested, returns (i, v). Otherwise returns the next 'v' only.
        Prints a progress bar (according to the fed estimated length)."""
        try:
            i, v = super().__next__()
        except StopIteration:
            time_passed_formatted = self._secs_to_days_hrs_mins_secs(self._time_stopper.get_total_runtime())
            print(self._progressbar_template.format(100, '00', '00', '00', '00', *time_passed_formatted), end='')
            raise

        # if not self._time_stopper.is_started():
        if i == 0:
            self._time_stopper.start_measure()
            print(self._progressbar_template.format(0, '--', '--', '--', '--', '--', '--', '--', '--'), end='')
        else:
            self._printing_time_counter += self._time_stopper.measure_time()
            # print a new time estimation & progress indicator every 3 seconds
            if self._printing_time_counter > self._update_frequency:
                # Main idea: give more weight to recent time estimations over old ones (e.g. maybe a program has
                # terminated/started, which affects our allocated CPU time.
                # time_estimation = time_per_itr * left_iters =
                # = (old_speed*0.95*old_iters + new_speed*recent_iters) / (0.95*old_iters + recent_iters) * left_iters
                # = ((old_time/old_iters)*0.95*old_iters + (recent_time/recent_iters)*recent_iters) / (0.95*old_iters + recent_iters) * left_iters
                # = (old_time*0.95 + recent_time) / (0.95*old_iters + recent_iters) * left_iters
                # Excplicitly:
                # recent_iters = i - self._left_time_estimation['num_of_iters']
                # new_time_estimation = ((0.95*self._left_time_estimation['passed_time'] + self._printing_time_counter) /
                #                        (0.95*self._left_time_estimation['num_of_iters'] + recent_iters) *
                #                        (self._estimated_length - i))
                new_time_estimation = ((0.8*self._left_time_estimation['passed_time'] + self._printing_time_counter) /
                                       (i - 0.2*self._left_time_estimation['num_of_iters']) *
                                       (self._estimated_length - i))

                total_time_passed = self._time_stopper.get_total_runtime()
                time_left_formatted = self._secs_to_days_hrs_mins_secs(new_time_estimation)
                time_passed_formatted = self._secs_to_days_hrs_mins_secs(total_time_passed)
                print(self._progressbar_template.format(int(100 * i / self._estimated_length), *time_left_formatted,
                                                        *time_passed_formatted), end='')
                sys.stdout.flush
                self._left_time_estimation['left_time'] = new_time_estimation
                self._left_time_estimation['num_of_iters'] = i
                self._left_time_estimation['passed_time'] = total_time_passed
                self._printing_time_counter = 0
        if self._should_enum:
            return i, v
        else:
            return v

    @staticmethod
    def _secs_to_days_hrs_mins_secs(seconds):
        # Converts seconds to int!
        minutes, seconds= divmod(round(seconds), 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        return days, hours, minutes, seconds


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)