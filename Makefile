CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude -I/tmp/int-c-cst/include -I/tmp/int-c-sheaf/include -I/tmp/int-c-hodge/include -I/tmp/int-c-ergodic/include
LDFLAGS = -lm

TESTS = test_cst_to_sheaf test_sheaf_to_hodge test_ergodic_to_cst \
        test_renormalization_cst test_west_african_sheaf test_full_pipeline

.PHONY: all test clean deps

all: $(TESTS)

test: $(TESTS)
	@echo ""
	@echo "═══════════════════════════════════════════════════"
	@echo "  integration-c: Eight libraries. One theorem."
	@echo "═══════════════════════════════════════════════════"
	@echo ""
	@pass=0; fail=0; total=0; \
	for t in $(TESTS); do \
		echo "── $$t ──"; \
		if ./$$t; then \
			pass=$$((pass + 1)); \
		else \
			fail=$$((fail + 1)); \
		fi; \
		total=$$((total + 1)); \
		echo ""; \
	done; \
	echo "═══════════════════════════════════════════════════"; \
	echo "  $$total tests: $$pass passed, $$fail failed"; \
	echo "═══════════════════════════════════════════════════"; \
	[ $$fail -eq 0 ]

test_%: tests/test_%.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TESTS)

deps:
	@echo "Headers expected at /tmp/int-c-{cst,sheaf,hodge,ergodic}/include/"
