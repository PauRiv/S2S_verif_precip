This folder contains the codes for the analysis conducted in "Assessment of S2S ensemble extreme precipitation forecasts over Europe", Rivoire et al.

The input data is the ECMWF's S2S precipitation hindcast data, cycle 47r2,  from 2001 to 2020.
It has 11 ensemble members, initialized twice a week with a lead time of 46 days.
The spatial domain here is Europe, [30degN ; 72degN] x [-15degE; 49.5degE].

I) get exceedences:

1. computate percentiles: "Get_percentile_seasonleadtime_dependant.R" and "Get_percentile_season_observation.R"

2. get exceedences over percentiles: "Get_exceedances_over_percentile.R" and "Get_exceedances_OBS.R"


II) Binary Loss Index:

3z. useful functions: "Binary_loss_function_afo_lead_member.R"
3a. compute the median member: "get_median_exceed_hind.R"
3b. compute the BLI: "Seasonal_BLF_for_median_memb.R" for the local analysis, "TempACC_seas_BLF_for_median_memb.R" for the temporal aggregation and "SpatiACC_seas_BLF_for_median_memb.R" for the spatial aggregation
3c. compute the confidence interval on the BLI: "Compute_seasonal_CI_BLS.R", "Compute_seas_CI_BLS_temporal_Accumu.R" and "Compute_seas_CI_BLS_spatial_Accumu.R"
3d. compute the last skilful day: "LastSkillfulDay_for_median_member.R"

III) Brier score:

4a. compute the Brier score: "Brier_score_local_seas.R", "Brier_score_tempAcc.R" and "Brier_score_spatiAcc.R"
4b. compute the climatological Brier score: "Brier_score_local_seas_climato.R", "Brier_score_tempAcc_climato.R" and "Brier_score_spatiAcc_climato.R"
4c. compute the last skillful day "BSS_last_days_0.R"

IV) Sensitivity analysis:
5a. get the nth member: "Select_nth_member.R"
