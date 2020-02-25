#%%
import csv
import os

import pandas as pd
import statsmodels.api as sm

def run_regressions(data, group_name, household_well_as):
    """For the distributed wells model and household wells model, return the regression results
    and append predicted urinary arsenic values to dataframe"""
    data, distributed_results = distributed_wells_regress(data)
    data, household_results = household_wells_regress(data, group_name, household_well_as)
    return distributed_results, household_results, data

#%% base case
def distributed_wells_regress(data):
    """For distributed wells model, return input data with added column for
    regression predicted values and return regression results."""
    # add a column of ones to the x-data to get a constant term in the model
    x = sm.add_constant(data.arsenic_ugl)

    model = sm.OLS(data.urine_as, x)
    results = model.fit()

    urine_as_pred = results.params[1]*data.arsenic_ugl + results.params[0]
    data['urine_as_pred_distributed'] = urine_as_pred
    return data, results

#%% add wells in family compound
def household_wells_regress(data, group_name, household_well_as):
    """For household wells model, return input data with added column for
    regression predicted values and return regression results."""
    # add a column of ones to the x-data to get a constant term in the model
    x = sm.add_constant(pd.concat([data.arsenic_ugl, data[household_well_as]], axis=1))

    model = sm.OLS(data.urine_as, x)
    results = model.fit()

    urine_as_pred = results.params[1]*data.arsenic_ugl + \
                    results.params[2]*data[household_well_as] + \
                    results.params[0]
    data['urine_as_pred_household'] = urine_as_pred
    with open(os.path.abspath('output_data/' + group_name +'_urine_as_pred_household.csv'), mode='w') as savefile:
        writer = csv.writer(savefile)
        data_for_csv = [data.arsenic_ugl, data['urine_as_pred_household']]
        writer.writerow(['arsenic_ugl','urine_as_pred_household'])
        writer.writerows(zip(*data_for_csv))
    return data, results