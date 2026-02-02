from ortools.sat.python import cp_model
import numpy as np
import time

# CONSTANTS #
# input your array of capacitor values here, they cannot excede 2 decimal points unless you
# change your scale factor
CAPACITORS = [  1, 
                4.7, 
                10, 
                47, 
                470, 
                220, 
                100, 
                330, 
                820, 
                1000]
# This is your +/- error range
ERROR = 3
# This is the net capacitor value target
TARGET = 13.81

class CapMatch():
  def __init__(self, capacitors, target, error):
    # model and solver
    self.model = cp_model.CpModel()
    self.solver = cp_model.CpSolver()
    
    # scale factors
    self.sf = 10 ** 2
    self.guard_sf = 10 ** 2
    self.total_sf = self.sf * self.guard_sf

    # initialize input vars
    self.raw_capacitors = capacitors
    self.raw_target = target
    self.raw_error = error
    self.capacitors, self.target, self.error = self.scale(capacitors, target, error)

    # limit to avoid int overflow
    self.limit = 10**12

    # penalties
    self.outside_range_p = 10**12
    self.direction_p = 10**6
    self.set_spread_w = 10**3 
    self.series_spread_w = 10 **3
    # self.accuracy_w = 10 ** 2
      
  def cap_setup(self):
    self.max_sum = sum(self.capacitors) 
    max_cap_value = max(self.capacitors)
    min_cap_value = min(self.capacitors)

    # initialize the indices
    self.idx_a = self.model.NewIntVar(0, len(self.capacitors) - 1, f'idx_a')
    self.idx_b = self.model.NewIntVar(0, len(self.capacitors) - 1, f'idx_b')
    self.idx_c = self.model.NewIntVar(0, len(self.capacitors) - 1, f'idx_c')
    self.idx_d = self.model.NewIntVar(0, len(self.capacitors) - 1, f'idx_d')

    # initialize the capacitors
    self.cap_a = self.model.NewIntVar(min_cap_value, max_cap_value, f'cap_a')
    self.cap_b = self.model.NewIntVar(min_cap_value, max_cap_value, f'cap_b')
    self.cap_c = self.model.NewIntVar(min_cap_value, max_cap_value, f'cap_c')
    self.cap_d = self.model.NewIntVar(min_cap_value, max_cap_value, f'cap_d')
    
    # assign cap to index
    self.model.AddElement(self.idx_a, self.capacitors, self.cap_a)
    self.model.AddElement(self.idx_b, self.capacitors, self.cap_b)
    self.model.AddElement(self.idx_c, self.capacitors, self.cap_c)
    self.model.AddElement(self.idx_d, self.capacitors, self.cap_d)
  
  def calc_high_res_equivalent_capacitance(self):
    self.total_sum =self.model.NewIntVar(1, self.max_sum, 'total_sum')
    self.model.Add(self.total_sum == (self.cap_a + self.cap_b + self.cap_c + self.cap_d))

    self.sum_ac = self.model.NewIntVar(1, self.limit, 'a_mult_b')
    self.sum_bd = self.model.NewIntVar(1, self.limit, 'c_mult_d')
    self.model.Add(self.sum_ac == self.cap_a + self.cap_c)
    self.model.Add(self.sum_bd == self.cap_b + self.cap_d)

    self.numerator = self.model.NewIntVar(0, self.limit, 'numerator')
    self.model.AddMultiplicationEquality(self.numerator, [self.sum_ac, self.sum_bd])

    self.scaled_numerator = self.model.NewIntVar(0, self.limit, 'scaled_numerator')
    self.model.Add(self.scaled_numerator == self.numerator * self.guard_sf)  

    self.rel_cap_high_res = self.model.NewIntVar(0, self.limit, 'rel_cap_high_res')
    self.model.AddDivisionEquality(self.rel_cap_high_res, self.scaled_numerator, self.total_sum)

    high_res_target = self.target * self.guard_sf
    self.high_res_error = self.error * self.guard_sf

    d = self.model.NewIntVar(-self.limit, self.limit, "d")
    self.model.Add(d == self.rel_cap_high_res - high_res_target)

    self.abs_d = self.model.NewIntVar(0, self.limit, 'abs_d')
    self.model.AddAbsEquality(self.abs_d, d)
  

  def calc_distance_penalties(self):
    # TODO: Normalize effective error
    self.effective_error = self.model.NewIntVar(0, self.limit, 'effective_err')

    # boolean for whether we are in or out of the error range
    self.is_outside = self.model.NewBoolVar('outside')

    # if the error_var is greater than the allowed error, then we are outside of the range
    self.model.Add(self.abs_d > self.high_res_error).OnlyEnforceIf(self.is_outside)
    self.model.Add(self.abs_d <= self.high_res_error).OnlyEnforceIf(self.is_outside.Not())

    # if inside, the distance from the target value is less important than minimizing the spread
    # if outside, the distance is more important
    self.model.Add(self.effective_error == self.abs_d * self.sf).OnlyEnforceIf(self.is_outside)
    self.model.Add(self.effective_error == self.abs_d).OnlyEnforceIf(self.is_outside.Not())

    # this makes sure that getting inside the error range is the most important goal, the penalty is huge outside, but 0 inside
    self.penalty_val = self.model.NewIntVar(0, self.limit, "penalty_val")
    self.model.Add(self.penalty_val == self.outside_range_p).OnlyEnforceIf(self.is_outside)
    self.model.Add(self.penalty_val == 0).OnlyEnforceIf(self.is_outside.Not())


  def normalize_set_spreads(self):
    # normalized spread: (a-c)/(a+c), spread is now a percentage of net capacitance
    self.spread_ac = self.model.NewIntVar(-self.limit, self.limit, 'spread_ac')
    self.model.Add(self.spread_ac == self.cap_a - self.cap_c)
    abs_spread_ac = self.model.NewIntVar(0, self.limit, 'spread_ac')
    self.model.AddAbsEquality(abs_spread_ac, self.spread_ac)
    scaled_spread_ac = self.model.NewIntVar(0, self.limit, 'scaled_spread_ac')
    self.model.Add(scaled_spread_ac == abs_spread_ac * self.guard_sf)
    self.norm_spread_ac = self.model.NewIntVar(0, self.limit, 'norm_spread_ac')
    self.model.AddDivisionEquality(self.norm_spread_ac, scaled_spread_ac, self.sum_ac)

    # normalized spread: (b-d)/(b+d), spread is now a percentage of net capacitance
    self.spread_bd = self.model.NewIntVar(-self.limit, self.limit, 'spread_bd')
    self.model.Add(self.spread_bd == self.cap_b - self.cap_d)
    abs_spread_bd = self.model.NewIntVar(0, self.limit, 'spread_bd')
    self.model.AddAbsEquality(abs_spread_bd, self.spread_bd)
    scaled_spread_bd = self.model.NewIntVar(0, self.limit, 'scaled_spread_bd')
    self.model.Add(scaled_spread_bd == abs_spread_bd * self.guard_sf)
    self.norm_spread_bd = self.model.NewIntVar(0, self.limit, 'norm_spread_bd')
    self.model.AddDivisionEquality(self.norm_spread_bd, scaled_spread_bd, self.sum_bd)
  
  def calc_parallel_bools(self):
    self.ac_par = self.model.NewBoolVar('ac_par')
    self.bd_par = self.model.NewBoolVar('bd_par')

    self.model.Add(self.cap_a == self.cap_c).OnlyEnforceIf(self.ac_par)
    self.model.Add(self.cap_a != self.cap_c).OnlyEnforceIf(self.ac_par.Not())
    self.model.Add(self.cap_b == self.cap_d).OnlyEnforceIf(self.bd_par)
    self.model.Add(self.cap_b != self.cap_d).OnlyEnforceIf(self.bd_par.Not())

  def calc_series_spread(self):
    self.series_spread = self.model.NewIntVar(-self.limit, self.limit, 'series_spread')
    self.model.Add(self.series_spread == self.sum_ac - self.sum_bd)

    abs_series_spread = self.model.NewIntVar(0, self.limit, 'series_spread')
    self.model.AddAbsEquality(abs_series_spread, self.series_spread)

    scaled_series_spread = self.model.NewIntVar(0, self.limit, 'scaled_series_spread')
    self.model.Add(scaled_series_spread == abs_series_spread * self.guard_sf)

    self.norm_series_spread = self.model.NewIntVar(0, self.limit, 'norm_series_spread')
    self.model.AddDivisionEquality(self.norm_series_spread, scaled_series_spread, self.total_sum)

    self.spread_is_positive = self.model.NewBoolVar('spread_is_positive')
    self.model.Add(self.series_spread > 0).OnlyEnforceIf(self.spread_is_positive)
    self.model.Add(self.series_spread <= 0).OnlyEnforceIf(self.spread_is_positive.Not())

    self.spread_is_negative = self.model.NewBoolVar('spread_is_negative')
    self.model.Add(self.series_spread < 0).OnlyEnforceIf(self.spread_is_negative)
    self.model.Add(self.series_spread >= 0).OnlyEnforceIf(self.spread_is_negative.Not())
         
  def calc_direction_penalty(self):
    is_violated = self.model.NewBoolVar('direction_violation')
    # violation if ac is parallel AND bd is mismatched AND ac sum > ac sum
    v1 = self.model.NewBoolVar('v1')
    self.model.AddBoolAnd([self.ac_par, self.bd_par.Not(), self.spread_is_positive]).OnlyEnforceIf(v1)
    self.model.AddBoolOr([self.ac_par.Not(), self.bd_par, self.spread_is_positive.Not()]).OnlyEnforceIf(v1.Not())

    # violation if bd is parallel AND ac is mismatched AND bd sum > ac sum
    v2 = self.model.NewBoolVar('v2')
    self.model.AddBoolAnd([self.ac_par.Not(), self.bd_par, self.spread_is_negative]).OnlyEnforceIf(v2)
    self.model.AddBoolOr([self.ac_par, self.bd_par.Not(), self.spread_is_negative.Not()]).OnlyEnforceIf(v2.Not())

    # link v1/v2 to the main violation boolean
    self.model.AddBoolOr([v1, v2]).OnlyEnforceIf(is_violated)
    self.model.AddBoolAnd([v1.Not(), v2.Not()]).OnlyEnforceIf(is_violated.Not())

    self.direction_penalty = self.model.NewIntVar(0, self.direction_p, 'direction_penalty')
    self.model.Add(self.direction_penalty == self.direction_p).OnlyEnforceIf(is_violated)
    self.model.Add(self.direction_penalty == 0).OnlyEnforceIf(is_violated.Not())

  def minimize(self):
    self.model.Minimize(
    self.penalty_val +                                  # must be in error range 
    self.direction_penalty +                            # satisfy direction rule    
    (self.set_spread_w * self.norm_spread_ac) +         # minimize mismatch 
    (self.set_spread_w * self.norm_spread_bd) +         # minimize mismatch 
    (self.series_spread_w * self.norm_series_spread) +  # balance ac vs bd sum
    (self.effective_error))                             # closeness to target

  def solve(self, debug=False):
  
    # set up the model solver
    if debug:
      self.solver.parameters.log_search_progress = True
      self.solver.parameters.log_to_stdout = True
      print(self.model.Validate())

    self.status = self.solver.Solve(self.model)

    if self.status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        self.print_report()
    else:
        print("No solution found.")

  def scale(self, capacitors, target, error):
    caps = (np.array(capacitors) * self.sf).astype(int).tolist()
    target_val = int(target * self.sf)
    error = int(error * self.sf)

    return caps, target_val, error

  def print_report(self):
    s = self.solver 
    # scaling back to pF
    a, b, c, d = s.Value(self.cap_a)/self.sf, s.Value(self.cap_b)/self.sf, s.Value(self.cap_c)/self.sf, s.Value(self.cap_d)/self.sf
    actual = s.Value(self.rel_cap_high_res) / self.total_sf

    print("="*35)
    print(f"{'MATCH NETWORK REPORT':^35}")
    print("="*35)
    print(f"Slot A: {a:<4} pF | Slot C: {c:<4} pF")
    print(f"Slot B: {b:<4} pF | Slot D: {d:<4} pF")
    print("-" * 35)
    print(f"Target: {self.raw_target} pF")
    print(f"Actual: {actual:.2f} pF")
    print(f"Error:  {actual - self.raw_target:.2f} pF")
    print("-" * 35)
    
    # Balance & Rules
    print(f"AC Mismatch: {s.Value(self.norm_spread_ac)/10}%")
    print(f"BD Mismatch: {s.Value(self.norm_spread_bd)/10}%")
    
    print("="*35)

  def print_single_value(self, name):
    s = self.solver
    if not hasattr(self, name):
        print(f"Error: Attribute '{name}' does not exist in MatchNetwork.")
        return
    
    var_handle = getattr(self, name)

    if not hasattr(self, 'status'):
        print("Error: Solver has not been run yet.")
        return

    raw_value = self.solver.Value(var_handle)

    if 'cap_' in name or 'sum_' in name:
        scaled_value = raw_value / self.sf
        unit = "pF"
    elif 'rel_cap' in name or 'abs_d' in name:
        scaled_value = raw_value / self.total_sf
        unit = "pF"
    else:
        scaled_value = raw_value
        unit = "(Raw)"

    print(f"{name: <25} : {scaled_value:<10} {unit}")
  
  def build_network(self):
    self.cap_setup()
    self.calc_high_res_equivalent_capacitance()
    self.calc_distance_penalties()
    self.normalize_set_spreads()
    self.calc_parallel_bools()
    self.calc_series_spread()
    self.calc_direction_penalty()
    self.minimize()
  


if __name__ == "__main__":
  start_time = time.perf_counter()

  network = CapMatch(CAPACITORS, TARGET, ERROR)
  network.build_network()
  network.solve()

  end_time = time.perf_counter()

  elapsed_time = end_time - start_time
  print(f"Execution time: {elapsed_time} seconds")
  