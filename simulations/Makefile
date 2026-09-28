# Run from the repository root with: make -C simulations
# This generates inputs for the 24 explicit PDMS chain-length/PDI cases.

.PHONY: all generate-all build-generator check-configs

all: generate-all

build-generator:
	$(MAKE) -C .. Generator/oil_generator

check-configs:
	@set -- N*_PDI*/model.conf; \
	if [ "$$#" -ne 24 ] || [ ! -f "$$1" ]; then \
		echo "Expected 24 N*_PDI*/model.conf files; found $$#" >&2; \
		exit 1; \
	fi

generate-all: check-configs build-generator
	@set -eu; \
	for config in N*_PDI*/model.conf; do \
		echo "Generating $$config"; \
		../Generator/oil_generator --config "$$config"; \
	done; \
	echo "Generated all 24 PDMS cases."
