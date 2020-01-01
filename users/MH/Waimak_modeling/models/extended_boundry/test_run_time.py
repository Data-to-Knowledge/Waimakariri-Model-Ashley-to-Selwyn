"""
 Author: Matt Hanson
 Created: 01-Jan-20 4:55 PM
 """

from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import get_stocastic_set
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.modpath_wrapper import import_gns_model
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.mt3d_wrapper import setup_run_mt3d, get_sft_stress_period_data, get_default_mt3d_kwargs, get_ssm_stress_period_data
from time import time
import os
import pandas as pd

if __name__ == '__main__':
    model = get_stocastic_set()[0]
    outdir = r"C:\matt_modelling_unbackedup\testmodflow_mt3d"
    out = {}
    t = time()
    m = import_gns_model(model,'test_run_time',outdir,mt3d_link=True)
    m.write_input()
    m.run_model()
    out['units']='minutes'
    out['modflow_run'] = (time()-t)/60
    print(out)

    ftl_path=os.path.join(m.model_ws,'{}.ftl'.format(m.name))

    t = time()
    ssm = get_ssm_stress_period_data(wil_race_con=1, upper_n_bnd_flux_con=0, lower_n_bnd_flux_con=0,
                               well_stream_con=0, llrzf_con=0, ulrzf_con=0, s_race_con=1, chb_con=0)

    sft = get_sft_stress_period_data(eyre=0, waimak=1, ash_gorge=1, cust=0, cust_biwash=1, ash_tribs=0,
                               glen_tui=None, garry=None, bullock=None, okuku=None, makerikeri=None, eyre_mar=0)

    setup_run_mt3d(ftl_path=ftl_path,mt3d_name='waimat_emma',mt3d_ws=outdir,ssm_crch=0,ssm_stress_period_data=ssm, sft_spd=sft)
    out['emma_run'] = (time()-t)/60
    print(out)


    t = time()
    sft = get_sft_stress_period_data()
    ssm = get_ssm_stress_period_data()
    setup_run_mt3d(ftl_path=ftl_path,mt3d_name='med_n',mt3d_ws=outdir,ssm_crch=0,ssm_stress_period_data=ssm, sft_spd=sft)
    out['med_n_run'] = (time()-t)/60
    print(out)

    out = pd.DataFrame(out)
    out.to_csv(os.path.join(outdir,'times.csv'))

