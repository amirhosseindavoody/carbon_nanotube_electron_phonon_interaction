
FC=gfortran
FCFLAGS = -O3 -ffree-line-length-none -fbounds-check
# FCFLAGS = -g
# FCFLAGS += -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wc-binding-type -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Wreal-q-constant -Wunused -Wconversion-extra -Wimplicit-interface -Wintrinsic-shadow -Wsurprising -Wfunction-elimination -fbounds-check -ffree-line-length-none -fall-intrinsics

SRCDIR = ./src
OBJDIR = ./obj

cnt_class.o: constants_mod.o math_functions_mod.o
cnt_electron_mod.o: cnt_class.o constants_mod.o graphene_mod.o math_functions_mod.o write_log_mod.o
cnt_geometry_mod.o: cnt_class.o constants_mod.o math_functions_mod.o write_log_mod.o
cnt_phonon_mod.o: cnt_class.o graphene_mod.o write_log_mod.o
cnt_scattering_electron_phonon_mod.o: cnt_class.o cnt_phonon_mod.o constants_mod.o graphene_mod.o math_functions_mod.o sim_properties_mod.o write_log_mod.o
cnt_scattering_exciton_phonon_mod.o: cnt_class.o cnt_phonon_mod.o constants_mod.o graphene_mod.o math_functions_mod.o sim_properties_mod.o write_log_mod.o
first_order_coulomb_transition_mod.o: constants_mod.o cnt_class.o geometric_matrix_element_mod.o kspace_matrix_element_mod.o math_functions_mod.o partition_function_mod.o rotate_shift_mod.o sim_properties_mod.o transition_points_mod.o write_log_mod.o
geometric_matrix_element_mod.o: cnt_class.o constants_mod.o
graphene_mod.o: constants_mod.o math_functions_mod.o
input_cnt_mod.o: cnt_class.o constants_mod.o write_log_mod.o
kspace_matrix_element_mod.o: cnt_class.o constants_mod.o write_log_mod.o
main.o: cnt_class.o cnt_electron_mod.o cnt_geometry_mod.o cnt_phonon_mod.o cnt_scattering_electron_phonon_mod.o cnt_scattering_exciton_phonon_mod.o first_order_coulomb_transition_mod.o input_cnt_mod.o partition_function_mod.o second_order_coulomb_phonon_transition_mod.o sim_properties_mod.o write_log_mod.o
partition_function_mod.o: cnt_class.o constants_mod.o sim_properties_mod.o
rotate_shift_mod.o: cnt_class.o math_functions_mod.o constants_mod.o
second_order_coulomb_phonon_transition_mod.o: constants_mod.o cnt_class.o sim_properties_mod.o write_log_mod.o
sim_properties_mod.o: cnt_class.o constants_mod.o write_log_mod.o
transition_points_mod.o: cnt_class.o constants_mod.o sim_properties_mod.o write_log_mod.o

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
