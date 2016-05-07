
FC=gfortran
# FCFLAGS = -O3 -ffree-line-length-none
FCFLAGS = -g
# FCFLAGS += -fcheck=all -Wall -Wconversion-extra -Wtabs -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination -fbounds-check -ffree-line-length-none -fall-intrinsics
FCFLAGS += -fcheck=all  -Waliasing -Wampersand -Wconversion -Wsurprising -Wc-binding-type -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Wreal-q-constant -Wunused -Wconversion-extra -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

# a2a_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# a2ep_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# a2em_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
cnt_class.o: constants_mod.o math_functions_mod.o
cnt_electron_mod.o: cnt_class.o constants_mod.o graphene_mod.o math_functions_mod.o write_log_mod.o
cnt_geometry_mod.o: cnt_class.o constants_mod.o math_functions_mod.o write_log_mod.o
cnt_phonon_mod.o: cnt_class.o graphene_mod.o write_log_mod.o
cnt_scattering_electron_phonon_mod.o: cnt_class.o cnt_phonon_mod.o constants_mod.o graphene_mod.o math_functions_mod.o sim_properties_mod.o write_log_mod.o
comparams.o: cnt_class.o
# em2a_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# em2ep_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# em2em_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# ep2a_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# ep2ep_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# ep2em_kspace_matrix_element_mod.o: comparams.o constants_mod.o transition_points_mod.o write_log_mod.o
# geometric_matrix_element_mod.o: comparams.o constants_mod.o
graphene_mod.o: constants_mod.o math_functions_mod.o
input_cnt_mod.o: cnt_class.o constants_mod.o write_log_mod.o
# main.o: cnt_class.o comparams.o input_cnt_mod.o occupation_mod.o parse_input_file_mod.o prepareForster_module.o sim_properties_mod.o transition_table_mod.o write_log_mod.o
main.o: cnt_class.o cnt_electron_mod.o cnt_geometry_mod.o cnt_phonon_mod.o cnt_scattering_electron_phonon_mod.o comparams.o input_cnt_mod.o sim_properties_mod.o write_log_mod.o
# occupation_mod.o: comparams.o cnt_class.o constants_mod.o
# parse_input_file_mod.o: comparams.o constants_mod.o transition_table_mod.o write_log_mod.o
# prepareForster_module.o: cnt_class.o comparams.o constants_mod.o
# rotate_shift_mod.o: cnt_class.o math_functions_mod.o constants_mod.o
sim_properties_mod.o: constants_mod.o write_log_mod.o
# transition_points_mod.o: cnt_class.o comparams.o math_functions_mod.o constants_mod.o write_log_mod.o
# transition_table_mod.o: a2a_kspace_matrix_element_mod.o a2ep_kspace_matrix_element_mod.o a2em_kspace_matrix_element_mod.o cnt_class.o comparams.o ep2a_kspace_matrix_element_mod.o ep2em_kspace_matrix_element_mod.o ep2ep_kspace_matrix_element_mod.o em2a_kspace_matrix_element_mod.o em2em_kspace_matrix_element_mod.o em2ep_kspace_matrix_element_mod.o geometric_matrix_element_mod.o constants_mod.o prepareForster_module.o rotate_shift_mod.o transition_points_mod.o write_log_mod.o

main: main.o
	$(FC) $(FCFLAGS) -o $@.exe $(wildcard $(OBJDIR)/*.o) -llapack -lblas
	@rm -f *.o *.mod
	@rm -rf $(OBJDIR)

%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) -c $(FCFLAGS) $< -J$(OBJDIR)
	@mv -f $@ $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.mod *.exe
	@rm -rf $(OBJDIR)
