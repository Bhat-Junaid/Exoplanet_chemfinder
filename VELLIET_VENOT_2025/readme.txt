R. Veillet     November 2024

Veillet et al., 2025, A&A

This network contains 1689 chemical reactions (1609 reversible and 80 irreversible). The reversal of reaction rates is based on the use of equilibrium constants.
These reactions involve 227 neutral species, that are listed in composes.dat.
(Same as Chemkin format, but without Argon.)

November 2025 Edit :
We have been informed of a typo in the Glarborg et al. 2014 kinetic network for the reaction NH2 + NH2 = N2H2 + H2, where n is 1.62 instead of 1.02 (Klippenstein et al. 2009).
This error was passed on to the present network, and was not fixed in this version to keep it identical to the one used in the paper.
This will be fixed in the next release, but it should be kept in mind when using the current version in simulations.
We also added a warning detailed below for the SO2 + H = SO + OH reaction that has a different formalism from the other Troe reactions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Combination reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

=================================================

The format of the reactions_1.dat file (Reversible combination reactions - Troe formalism, combinaison_k0_kinf_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    10(1x,e10.3)                                                               4(1x,e10.3)                    17(1x,e10.3)
         Reactants    Products     Alpha_0 Beta_0 Gamma_0 F_0 g_0 Alpha_inf Beta_inf Gamma_inf F_inf g_inf    A_Troe B_Troe C_Troe D_Troe    Efficiencies

Alpha_0, Beta_0, Gamma_0 are the parameters to compute the rate coefficients in the low pressure limits k_0. Alpha_inf, Beta_inf, Gamma_inf are the parameters to compute the rate coefficients in the low pressure limits k_inf. The formula adopted follow the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F_0 is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g_0 is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

A_Troe, B_Troe, C_Troe, D_Troe are the four parameters needed to compute the fall-off Troe formulation. When these four parameters are set equal to 0, the Lindemann formulation is adopted to compute the fall-off (see http://kida.obs.u-bordeaux1.fr/help).

Enhanced third body efficiencies are specified for the following selected species :
CH3OCH3, He, H2S, H, SO2, CO, N2, O2, CH4, C2H6, CO2, C2H4, AR, H2O, H2, H2O2, C2H2

/!\ WARNING /!\

For reactions with two reactants and two products in this file (currently only SO2 + H = SO + OH), the rate constant expression is different.
These reactions (chemically activated reactions) follow another formalism that can be found at https://www.cantera.org/stable/reference/kinetics/rate-constants.html
Instead of the usual k = F*k_inf*P_r/(1 + P_r), for these reactions we have k = F*k_0/(1 + P_r).

This also results in a different asymptotic behavior as we have :
k -> k_0       when P -> 0
k -> k_inf/[M] when P -> +inf
instead of the usual :
k -> k_0*[M]   when P -> 0
k -> k_inf     when P -> +inf

This also means that the units of k_inf and k_0 must be different from the usual Troe formalism :
k_0 : cm³/molecule/s
k_inf : 1/s
instead of the usual :
k_0 : cm⁶/molecule²/s
k_inf : cm³/molecule/s

/!\ WARNING /!\

=================================================

The format of the reactions_2.dat file (Reversible combination reactions - low pressure limit, combinaison_k0_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    5(1x,e10.3)                       17(1x,e10.3)
         Reactants    Products     Alpha_0 Beta_0 Gamma_0 F_0 g_0    Efficiencies

Alpha_0, Beta_0, Gamma_0 are the parameters to compute the rate coefficients in the low pressure limits k_0. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F_0 is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g_0 is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

Enhanced third body efficiencies are specified for the following selected species :
CH3OCH3, He, H2S, H, SO2, CO, N2, O2, CH4, C2H6, CO2, C2H4, AR, H2O, H2, H2O2, C2H2

=================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Decomposition reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

=================================================

The format of the reactions_3.dat file (Irreversible decomposition reactions, decompo_irrev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    5(1x,e10.3)
         Reactants    Products     Alpha Beta Gamma F g

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

The format of the reactions_4.dat file (Irreversible decomposition reactions - coefficients with a pressure dependence, decompo_irrev_plog.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    2(1x,e10.3)
         Reactants    Products     P              Alpha          Beta           Gamma          F g
with N an integer corresponding to the number of pressures and their corresponding set of rate constant parameters, used to interpolate the rate constant pressure dependance.

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help
Several set of the parameters Alpha, Beta, Gamma are given, corresponding to different pressures. The corresponding pressures are given in the column ‘P’, in atm. To be more precise, the first Alpha, Beta, Gamma of each column correspond to the first value of P. 
To calculate the reaction rate in your model, the rate parameters are linearly interpolated with respect to the logarithm of the pressure in between two pressure.  The rate coefficients will take the lowest pressure value or the highest pressure value when the pressure is outside the specified range. 

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

The format of the reactions_5.dat file (Reversible decomposition reactions - Troe formalism, decompo_k0_kinf_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    10(1x,e10.3)                                                               4(1x,e10.3)                    17(1x,e10.3)
         Reactants    Products     Alpha_0 Beta_0 Gamma_0 F_0 g_0 Alpha_inf Beta_inf Gamma_inf F_inf g_inf    A_Troe B_Troe C_Troe D_Troe    Efficiencies

Alpha_0, Beta_0, Gamma_0 are the parameters to compute the rate coefficients in the low pressure limits k_0. Alpha_inf, Beta_inf, Gamma_inf are the parameters to compute the rate coefficients in the low pressure limits k_inf. The formula adopted follow the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F_0 is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g_0 is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

A_Troe, B_Troe, C_Troe, D_Troe are the four parameters needed to compute the fall-off Troe formulation. When these four parameters are set equal to 0, the Lindemann formulation is adopted to compute the fall-off (see http://kida.obs.u-bordeaux1.fr/help).

Enhanced third body efficiencies are specified for the following selected species :
CH3OCH3, He, H2S, H, SO2, CO, N2, O2, CH4, C2H6, CO2, C2H4, AR, H2O, H2, H2O2, C2H2

=================================================

The format of the reactions_6.dat file (Reversible decomposition reactions - low pressure limit, decompo_k0_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    5(1x,e10.3)                       17(1x,e10.3)
         Reactants    Products     Alpha_0 Beta_0 Gamma_0 F_0 g_0    Efficiencies

Alpha_0, Beta_0, Gamma_0 are the parameters to compute the rate coefficients in the low pressure limits k_0. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F_0 is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g_0 is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

Enhanced third body efficiencies are specified for the following selected species :
CH3OCH3, He, H2S, H, SO2, CO, N2, O2, CH4, C2H6, CO2, C2H4, AR, H2O, H2, H2O2, C2H2

=================================================

The format of the reactions_7.dat file (Reversible decomposition reactions, decompo_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    5(1x,e10.3)
         Reactants    Products     Alpha Beta Gamma F g

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

The format of the reactions_8.dat file (Reversible decomposition reactions - coefficients with a pressure dependence, decompo_rev_plog.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    2(1x,e10.3)
         Reactants    Products     P              Alpha          Beta           Gamma          F g
with N an integer corresponding to the number of pressures and their corresponding set of rate constant parameters, used to interpolate the rate constant pressure dependance.

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help
Several set of the parameters Alpha, Beta, Gamme are given, corresponding to different pressures. The corresponding pressures are given in the column ‘P’, in atm. To be more precise, the first Alpha, Beta, Gamma of each column correspond to the first value of P. 
To calculate the reaction rate in your model, the rate parameters are linearly interpolated with respect to the logarithm of the pressure in between two pressure. The rate coefficients will take the lowest pressure value or the highest pressure value when the pressure is outside the specified range.

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Desexcitation reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

=================================================

The format of the reactions_9.dat file (Reversible third-body assisted desexcitation reactions, desexcitation_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    5(1x,e10.3)                       17(1x,e10.3)
         Reactants    Products     Alpha_0 Beta_0 Gamma_0 F_0 g_0    Efficiencies

Alpha_0, Beta_0, Gamma_0 are the parameters to compute the rate coefficients in the low pressure limits k_0. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F_0 is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g_0 is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

Enhanced third body efficiencies are specified for the following selected species :
CH3OCH3, He, H2S, H, SO2, CO, N2, O2, CH4, C2H6, CO2, C2H4, AR, H2O, H2, H2O2, C2H2

=================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bimolecular reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

=================================================

The format of the reactions_10.dat file (Irreversible bimolecular reactions, reaction_2_Corps_irrev.dat in FRECKLL format) is the following:

Format : 5(1x,a10) 5(1x,a10) 5(1x,e10.3)
         Reactants   Products  Alpha  Beta  Gamma  F g

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

The format of the reactions_11.dat file (Irreversible bimolecular reactions - coefficients with a pressure dependence, reaction_2_Corps_irrev_plog.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    2(1x,e10.3)
         Reactants    Products     P              Alpha          Beta           Gamma          F g
with N an integer corresponding to the number of pressures and their corresponding set of rate constant parameters, used to interpolate the rate constant pressure dependance.

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help
Several set of the parameters Alpha, Beta, Gamme are given, corresponding to different pressures. The corresponding pressures are given in the column ‘P’, in atm. To be more precise, the first Alpha, Beta, Gamma of each column correspond to the first value of P. 
To calculate the reaction rate in your model, the rate parameters are linearly interpolated with respect to the logarithm of the pressure in between two pressure.  The rate coefficients will take the lowest pressure value or the highest pressure value when the pressure is outside the specified range. 

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

The format of the reactions_12.dat file (Reversible bimolecular reactions, reaction_2_Corps_rev.dat in FRECKLL format) is the following:

Format : 5(1x,a10) 5(1x,a10) 5(1x,e10.3)
Reactants   Products  Alpha  Beta  Gamma  F g

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

The format of the reactions_13.dat file (Reversible bimolecular reactions - coefficients with a pressure dependence, reaction_2_Corps_rev_plog.dat in FRECKLL format) is the following:

Format : 5(1x,a10)    5(1x,a10)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    N(1x,e10.3)    2(1x,e10.3)
         Reactants    Products     P              Alpha          Beta           Gamma          F g
with N an integer corresponding to the number of pressures and their corresponding set of rate constant parameters, used to interpolate the rate constant pressure dependance.

Alpha, Beta, Gamma are the parameters to compute the rate coefficients. The formula adopted follows the Kooij formalism.
See http:// kida.obs.u-bordeaux1.fr/help
Several set of the parameters Alpha, Beta, Gamme are given, corresponding to different pressures. The corresponding pressures are given in the column ‘P’, in atm. To be more precise, the first Alpha, Beta, Gamma of each column correspond to the first value of P. 
To calculate the reaction rate in your model, the rate parameters are linearly interpolated with respect to the logarithm of the pressure in between two pressure.  The rate coefficients will take the lowest pressure value or the highest pressure value when the pressure is outside the specified range. 

F is the uncertainty factor on the rate coefficient (see http://kida.obs.u-bordeaux1.fr/help)
g is the temperature dependence of this uncertainty factor (see http://kida.obs.u-bordeaux1.fr/help)
Type of uncertainty is lognormal (logn)

=================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
