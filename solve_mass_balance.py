import csv
import os

import sympy as sym

def calculate_parameters(distributed_results, household_results, ext_params, group_name):
    """Calculate parameters for distributed well model and household well model
    and save to csv files."""
    solve_params_distributed(distributed_results, group_name, ext_params)
    solve_params_household(household_results, group_name, ext_params)

# adapted function from Jason
def propagate_uncertainty(expr, values=None, uncertainties=None, d=sym.symbols('d', cls=sym.Function)):
    """
    propagate_uncertainty takes a sympy expression, and optionally dicts mapping variables
    to their values and uncertainties, and returns the propagated uncertainty

    >>> x, y = sym.symbols('x y')
    >>> propagate_uncertainty(x * y)
    (x*y, sqrt(x**2*d(y)**2 + y**2*d(x)**2))
    >>> propagate_uncertainty(x * y, {x:2, y:3})
    (6, sqrt(9*d(x)**2 + 4*d(y)**2))
    >>> propagate_uncertainty(x * y, {x:20, y:30}, {x:2, y:3})
    (600, 60*sqrt(2))
    """
    if values is None: values = {}
    if uncertainties is None: uncertainties = {}
    variables = list(expr.free_symbols)
    dv = (lambda v: sym.symbols(f'd{v.name}'))
    # check that we aren't reusing variable names:
    for v in variables:
        if dv(v) in variables:
            raise ValueError(f'propagate_uncertainty cannot handle having both a variable named {(repr(str(v)))} and a variable named {repr(str(dv(v)))}')
    # calculate symbolic error, as per https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Simplification
    sym_error = sym.sqrt(sum(sym.diff(expr, v)**2 * dv(v)**2 for v in variables))
    # substitute in the uncertainties and the values
    # order of substitution does not matter because we use variables rather than
    #  functions for the differentials, and because we fail if any of the differential
    #  variable names are already used as variable names
    num_error = sym_error.subs({dv(k):v for (k, v) in uncertainties.items()}).subs(values)
    # replace the remaining uncertainties with derivative functions, rather than custom-named variables
    # this allows the user to specify the form the differentials take
    #  without resulting in us substituting in arguments to the differential function
    #  (that is, we don't want to return 9*d(2)**2 + 4*d(3)**2 for propagate_uncertainty(x * y, {x:2, y:3})
    final_error = num_error.subs({dv(k):d(k) for k in variables})
    return expr.subs(values), final_error

def split_uncertainties(parameters_with_uncertainties):
    """Takes a dict mapping variable to tuples of (value, uncertainty).
    Splits it into separate dicts for values and uncertainties and turns variables into sympy Symbols."""
    params = {sym.Symbol(k):v[0] for k, v in parameters_with_uncertainties.items()}
    uncertainties = {sym.Symbol(k):v[1] for k, v in parameters_with_uncertainties.items()}
    return params, uncertainties

def apply_formatting(arg):
    """Applies formatting to parameters for display."""
    if isinstance(arg, dict):
        return {k:apply_formatting(v) for k, v in arg.items()}
    elif isinstance(arg, tuple) and len(arg) == 2:
        val, uncertainty = arg
        return f'{float(val):.2f}+-{float(uncertainty):.2f}'
    else:
        return round(float(arg), 2)

def solve_params_distributed(model, group_name, ext_params):
    """Solve for parameters and uncertainties for distributed well model and save to csv."""
    ff, fc, md, mb, Mf, Q, avgAs, slope, intercept = sym.symbols('ff fc md mb Mf Q avgAs slope intercept')
    # get values and uncertainties for external parameters
    values, uncertainties = split_uncertainties(ext_params)
    # get values and uncertainties for slope and intercept
    values[slope] = model.params[1]
    values[intercept] = model.params[0]
    uncertainties[slope] = model.bse[1]
    uncertainties[intercept] = model.bse[0]

    # solve for fp, fu, and fo
    fp = (slope*(1-ff-fc)*avgAs+Mf/Q)/(slope*avgAs+intercept)
    fu = (1-md-mb)*(fp/slope)
    fo = 1 - fp - ff - fc
    frac_primary_well = fp/(fp+fo)
    frac_other_well = fo/(fp+fo)

    vfp, dfp = propagate_uncertainty(fp, values, uncertainties)
    vfu, dfu = propagate_uncertainty(fu, values, uncertainties)
    vfo, dfo = propagate_uncertainty(fo, values, uncertainties)
    vfrac_primary_well, dfrac_primary_well = propagate_uncertainty(frac_primary_well, values, uncertainties)
    vfrac_other_well, dfrac_other_well = propagate_uncertainty(frac_other_well, values, uncertainties)

    solutions = {'nobs':model.nobs, 'r2':model.rsquared, 'r2_adj':model.rsquared_adj,
                'intercept':(model.params[0], model.bse[0]), 'slope':(model.params[1], model.bse[1]),
                'fu':(vfu, dfu), 'fp':(vfp, dfp), 'fo':(vfo, dfo),
                'frac_primary_well':(vfrac_primary_well, dfrac_primary_well),
                'frac_other_well':(vfrac_other_well, dfrac_other_well)}
    # also include the input parameters so they can be referenced alongside the output parmeters
    solutions.update(ext_params)
    with open(os.path.abspath('output_data/' + group_name + '_distributed_solved.csv'), "w") as savefile:
        writer = csv.writer(savefile)
        for row in apply_formatting(solutions).items():
            writer.writerow(row)

def solve_params_household(model, group_name, ext_params):
    """Solve for parameters and uncertainties for distributed well model and save to csv."""
    ff, fc, md, mb, Mf, Q, avgAs, slope_primary, slope_household, intercept = \
        sym.symbols('ff fc md mb Mf Q avgAs slope_primary slope_household intercept')
    values, uncertainties = split_uncertainties(ext_params)

    # get values and uncertainties for slope and intercept
    values[slope_primary] = model.params[1]
    values[slope_household] = model.params[2]
    values[intercept] = model.params[0]
    uncertainties[slope_primary] = model.bse[1]
    uncertainties[slope_household] = model.bse[2]
    uncertainties[intercept] = model.bse[0]

    # solve for fp, fu, fo, fh
    fu = (1 - md - mb)*(1 - ff - fc + Mf/Q/avgAs)/(slope_primary + slope_household + intercept/avgAs)
    fp = slope_primary*(1 - ff - fc + Mf/Q/avgAs)/(slope_primary + slope_household + intercept/avgAs)
    fo = intercept*(1 - ff - fc + Mf/Q/avgAs)/(slope_primary*avgAs + slope_household*avgAs + intercept
        ) - Mf/Q/avgAs
    fh = 1 - fp - fo - ff - fc
    frac_primary_well = fp/(fp+fo+fh)
    frac_other_well = fo/(fp+fo+fh)
    frac_household_well = fh/(fp+fo+fh)

    vfp, dfp = propagate_uncertainty(fp, values, uncertainties)
    vfu, dfu = propagate_uncertainty(fu, values, uncertainties)
    vfo, dfo = propagate_uncertainty(fo, values, uncertainties)
    vfh, dfh = propagate_uncertainty(fh, values, uncertainties)
    vfrac_primary_well, dfrac_primary_well = propagate_uncertainty(frac_primary_well, values, uncertainties)
    vfrac_other_well, dfrac_other_well = propagate_uncertainty(frac_other_well, values, uncertainties)
    vfrac_household_well, dfrac_household_well = propagate_uncertainty(frac_household_well, values, uncertainties)

    solutions =  {'nobs':model.nobs, 'r2':model.rsquared, 'r2_adj':model.rsquared_adj,
                'intercept':(model.params[0], model.bse[0]),
                'slope_primary':(model.params[1], model.bse[1]),
                'slope_household':(model.params[2], model.bse[2]),
                'fu':(vfu, dfu), 'fp':(vfp, dfp), 'fh':(vfh, dfh), 'fo':(vfo, dfo),
                'frac_primary_well':(vfrac_primary_well, dfrac_primary_well),
                'frac_household_well':(vfrac_household_well, dfrac_household_well),
                'frac_other_well':(vfrac_other_well, dfrac_other_well)
                }
    solutions.update(ext_params)
    with open(os.path.abspath('output_data/' + group_name + '_household_solved.csv'), "w") as savefile:
        writer = csv.writer(savefile)
        for row in apply_formatting(solutions).items():
            writer.writerow(row)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
