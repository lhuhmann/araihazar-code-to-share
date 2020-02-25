#%% imports
import os
import pandas as pd

from regressions import run_regressions
from plots import make_plots
from solve_mass_balance import calculate_parameters

#%%
def make_subset(data, group_name): # pylint: disable=too-many-return-statements
    """Given a group name and the full data set, return the subset of the data that corresponds
    to the group name."""
    if group_name == 'did_not_know':
        # keep only the rows where people did not know their well As
        # when their urine As was measured
        return data[data['knew_well_as'] == False]
    if group_name == 'may_have_known':
        # keep only the rows where people may have known their well As
        # when their urine As was measured
        return data[data['knew_well_as'] == True]
    if group_name == 'women':
        return data[data['sex'] == 'female']
    if group_name == 'women_did_not_know':
        return data[(data['knew_well_as'] == False) & (data['sex'] == 'female')]
    if group_name == 'men':
        return data[data['sex'] == 'male']
    if group_name == 'men_did_not_know':
        return data[(data['knew_well_as'] == False) & (data['sex'] == 'male')]
    if group_name == 'all':
        return data
    assert False, 'Group does not exist'
    return False

def run_all():
    # parameters that go into the mass balance equation and their uncertainties
    parameters_with_uncertainties = {'ff':(0.2, 0.1), 'fc':(0.12, 0.06), 'md':(0.06, 0.03),
                                     'mb':(0, 0), 'Mf':(64, 4), 'Q':(3, 1), 'avgAs':(95.2, 1.4)}

    # which column to use for neighbor well arsenic
    household_well_as = "other_as_50m"
    # number of bins
    numbins = 15

    data = pd.read_csv(os.path.abspath("data_for_regressions.csv"))
    group_name = 'all'

    # get the correct subset of the data
    data = make_subset(data, group_name)
    # run regressions
    simple_results, two_slope_results, both_data = run_regressions(data, group_name, household_well_as)
    # calculate parameter values
    calculate_parameters(simple_results, two_slope_results, parameters_with_uncertainties, group_name)
    # plot results
    make_plots(simple_results, two_slope_results, both_data, group_name, numbins, household_well_as)

if __name__ == "__main__":
    run_all()
