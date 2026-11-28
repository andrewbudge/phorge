CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Werror -D_POSIX_C_SOURCE=200809L
LDLIBS = -lm -lminizip -lz
RUSTC ?= rustc
RUSTFLAGS ?= -O --edition 2024
VALGRIND ?= valgrind
VALGRIND_FLAGS ?= --leak-check=full --show-leak-kinds=definite,indirect --errors-for-leak-kinds=definite,indirect --error-exitcode=99
PHORGE_RUST ?= reference/phorge
CLIPKIT ?= /home/andrew-budge/miniconda3/envs/clipkit/bin/clipkit

.PHONY: all install release rust-reference golden-capture test test-gzip test-datasets test-rip valgrind bench sizes clipkit bench-clipkit bench-simd clean

PREFIX ?= /usr/local
BUILD := build
SRC := src
TESTS := tests
CSRC := $(wildcard $(SRC)/*.c)
HDRS := $(wildcard $(SRC)/*.h)
SANFLAGS := -std=c11 -g -O1 -Wall -Wextra -Werror -D_POSIX_C_SOURCE=200809L -fsanitize=address,undefined

all: $(BUILD)/phorge $(BUILD)/test_c

install: $(BUILD)/phorge
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(BUILD)/phorge $(DESTDIR)$(PREFIX)/bin/phorge

release:
	bash scripts/release.sh

golden-capture:
	@python3 bench/gen_fasta.py > /dev/null
	@python3 bench/gen_extract.py > /dev/null
	python3 bench/diff_test.py --capture --rust $(PHORGE_RUST)
	python3 bench/gzip_test.py --capture --rust $(PHORGE_RUST)

rust-reference:
	cargo build --release --manifest-path rust-legacy/Cargo.toml
	mkdir -p reference
	cp rust-legacy/target/release/phorge reference/phorge

$(BUILD):
	mkdir -p $@

$(BUILD)/phorge: $(CSRC) $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) $(CSRC) $(LDLIBS) -o $@

$(BUILD)/test_c: $(TESTS)/test_parsing.c $(SRC)/phorge_c.c $(SRC)/phorge_c.h $(HDRS) | $(BUILD)
	$(CC) $(SANFLAGS) -I$(SRC) $(TESTS)/test_parsing.c $(SRC)/phorge_c.c $(LDLIBS) -o $@

$(BUILD)/test_c_plain: $(TESTS)/test_parsing.c $(SRC)/phorge_c.c $(SRC)/phorge_c.h $(HDRS) | $(BUILD)
	$(CC) -std=c11 -g -O1 -Wall -Wextra -Werror -D_POSIX_C_SOURCE=200809L -I$(SRC) \
		$(TESTS)/test_parsing.c $(SRC)/phorge_c.c $(LDLIBS) -o $@

$(BUILD)/phorge-asan: $(CSRC) $(HDRS) | $(BUILD)
	$(CC) $(SANFLAGS) -I$(SRC) $(CSRC) $(LDLIBS) -o $@

test: all
	@python3 bench/gen_fasta.py > /dev/null
	@python3 bench/gen_extract.py > /dev/null
	ASAN_OPTIONS=detect_leaks=1 ./$(BUILD)/test_c
	python3 bench/diff_test.py
	$(MAKE) $(BUILD)/phorge-asan
	ASAN_OPTIONS=detect_leaks=1 python3 bench/diff_test.py --c ./$(BUILD)/phorge-asan
	ASAN_OPTIONS=detect_leaks=1 python3 bench/diff_test.py --c ./$(BUILD)/phorge-asan --only-extra
	python3 bench/gzip_test.py --c ./$(BUILD)/phorge
	ASAN_OPTIONS=detect_leaks=1 python3 bench/gzip_test.py --c ./$(BUILD)/phorge-asan
	$(MAKE) test-datasets
	$(MAKE) test-rip
	@if [ -x "$(CLIPKIT)" ]; then python3 bench/clipkit_diff.py --clipkit "$(CLIPKIT)" --rust $(PHORGE_RUST) --c ./$(BUILD)/phorge; else echo "skip: clipkit not found at $(CLIPKIT)"; fi

test-gzip: all
	python3 bench/gzip_test.py --c ./$(BUILD)/phorge

test-datasets: $(BUILD)/phorge
	python3 bench/datasets_test.py --c ./$(BUILD)/phorge

test-rip: $(BUILD)/phorge
	python3 bench/seqripper_test.py --c ./$(BUILD)/phorge
	python3 tests/rip_aggregate_test.py --c ./$(BUILD)/phorge

valgrind: $(BUILD)/phorge $(BUILD)/test_c_plain test-datasets
	VALGRIND="$(VALGRIND)" VALGRIND_FLAGS="$(VALGRIND_FLAGS)" \
		bash bench/valgrind_test.sh ./$(BUILD)/phorge ./$(BUILD)/test_c_plain

clipkit: all
	python3 bench/clipkit_diff.py --clipkit "$(CLIPKIT)" --rust $(PHORGE_RUST) --c ./$(BUILD)/phorge

bench-clipkit: all
	bash bench/run_clipkit_bench.sh

bench-simd: all
	bash bench/run_simd_bench.sh

bench: all
	python3 bench/gen_fasta.py
	bash bench/run_bench.sh

sizes: all
	bash bench/sizes.sh

clean:
	rm -rf $(BUILD)
	rm -rf bench/data/_out
