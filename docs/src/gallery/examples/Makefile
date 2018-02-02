EXAMPLE_SYSTEMS = \
henon_map \
standard_map \
bakers_map \
arnold_cat_map \
linz_sprott_99 \
van_der_pol \
beer_95 \
# lorenz_63

EXAMPLE_SCRIPTS = $(patsubst %,%.jl,$(EXAMPLE_SYSTEMS))

all: $(EXAMPLE_SCRIPTS)

$(EXAMPLE_SCRIPTS): %.jl: lorenz_63.jl
	sed -E 's/(LyapunovExponents)\.lorenz_63/\1.$*/' $^ > $@
