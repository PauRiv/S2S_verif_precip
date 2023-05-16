This folder contains the codes for the analysis conducted in "Assessment of S2S ensemble extreme precipitation forecast skill over Europe", Rivoire et al. 2022, submitted to HNESS

The input data is the ECMWF's S2S precipitation hindcast data, cycle 47r2,  from 2001 to 2020.
It has 11 ensemble members, initialized twice a week with a lead time of 46 days.
The spatial domain here is Europe, [30degN ; 72degN] x [-15degE; 49.5degE].

I) get percentiles and exceedences: R codes

1_Get_binary_extremes_OBS.R
2_Get_binary_extremes_HIND.R


II) Binary Loss Index: R codes

Binary_loss_index_functions.R: usefull functions

1_Get_median_member_hind.R
2_BLI_local_seasonal.R
3_BLI_temporally_aggreg_seasonal.R
4_BLI_spatially_aggreg_seasonal.R
5_BLI_spatially_aggreg_NAOs.R

6_confidence_interval_seasonal_BLI.R
7_confidence_interval_seasonal_temp_accum_BLI.R
8_confidence_interval_seasonal_spa_accum_BLI.R
9_confidence_interval_seasonal_NAOs_spa_accum_BLI.R
10_Last_skillful_day_BLI.R


III) Brier score: R codes

1_Brier_score_local_seasonal.R
2_Brier_score_local_seasonal_climatology.R
3_Last_skillful_day_BrierScore.R

IV) Figures: python notebook
Plot_last_skillful_day.ipynb
