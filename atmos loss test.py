from Calculators.function_solvers import calculate_total_mass_loss
from astropy import units
from astropy.constants import R_sun, L_sun, sigma_sb, M_earth, R_earth
from Data_management import nasa_exoplanet_download
from Data_management.file_manager import save_catalog
print("test for TOI-5734 b")
loss = calculate_total_mass_loss(0.1, (9.1*M_earth).to(units.kg), (2.1*1.4*R_earth).to(units.m), (2.1*R_earth).to(units.m), 10**-0.74232*3.828*10**26*units.W, 0.001, 0.1*1e9*365*24*60*60*units.s, 10*1e9*365*24*60*60*units.s, 1.23, (0.05921*units.AU).to(units.m))
print(loss*(9.1*M_earth).to(units.kg))
loss2 = calculate_total_mass_loss(0.1, (9.1*M_earth).to(units.kg), (2.1*1.4*R_earth).to(units.m), (2.1*R_earth).to(units.m), 10**-0.74232*3.828*1e26*units.W, 0.001, (0.1*units.Gyr).to(units.s), (10*units.Gyr).to(units.s), 1.23, (0.05921*units.AU).to(units.m))
print(loss2*(9.1*M_earth).to(units.kg))
print(0.6*1e9*365*24*60*60*units.s)
print((0.6*units.gigayear).to(units.s))
print(loss.decompose())
print(loss2.decompose())

catalog = nasa_exoplanet_download.main()
loss = calculate_total_mass_loss(0.1, (catalog["pl_masse"]).to(units.kg), (catalog["pl_rade"]*1.4).to(units.m), (catalog["pl_rade"]).to(units.m), catalog["st_lum"]*L_sun, catalog["Lxuv/Lbol"], catalog["t_sat"].to(units.s), (10*units.Gyr).to(units.s), 1.23, catalog["pl_orbsmax"].to(units.m))
#print(loss)

t_end_gyr = (10*units.Gyr).to_value(units.Gyr)

colname = f"% lost at t={t_end_gyr:.1f} Gyr"
catalog[colname] = loss.decompose() * 100
catalog[colname].format = ".3f"
save_catalog(catalog, "AR")

print(loss.decompose())
