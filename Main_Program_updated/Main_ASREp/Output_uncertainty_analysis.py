# %% IMPORTANT!!! READ THIS ENTIRE FILE:
#
# In here the output wanted can be specified:
# "error variable" is mandatory to include in the write list AS THE LAST VARIABLE
# The rest of the variables can be selected from the list below:


# output_fields = ['max_exx', 'e_t_b', 'max_e_xy', 'max_e_t_mid', 'max_e_t', 'error variable']
output_fields = ['max_exx', 'e_t_b', 'max_e_xy', 'max_e_t_mid', 'max_e_t', 'max_eps_t_ssi', 'error variable']


"""  VARIABLES TO INCLUDE: MORE CAN BE ADDED IN "FEA_LTS" in the function "compute_tensile_strain_Jinyan"
     SO FAR ONLY VARIABLES WITH 1 number has been validated, but it should be possible to extract vectors
     Remember, the 'error variable' must be last in 'output_fields' and the way the string is made is THE WAY IT WILL 
     BE PASSED TO QUOFEM
        'max_exx'                   Maximum longitudinal normal strain
        'exx,t,b':                  tensile_strain_top due to bending,
        'exx,b,b':                  tensile_strain_bottom due to bending,
        'e_t_b':                    Maximum tensile strains due to bending
        'max_e_xy':                 Maximum true_shear_strain,
        'max_e_t_mid':              Maximum tensile_strain_midpoint of beam,
        'max_e_t':                  Maximum tensile stain over the entire beam
        'max_eps_t_ssi'             Highest LTS damage category over the beam
        'error variable'            Specifies if results should be removed from quoFEM due to an error

"""

