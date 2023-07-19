# QUEST script file for getting KL xml statistics out
# assume we're running this in the problem directory

cd(/Users/gbedonia/temphumid_cc/sboom_atm_uq_temphumid_dcc_kl5l3)
remote(false)
panel(xml) read_database(*) read_measurement(TEMP_profile.bundle.xml) \
            import_curve(TEMPMEAN_profile.txt) \
            read_measurementconfig(configtempstatextract) \
            write_data(/stats/KLTEMP) 
