

install: find_nonsynon_mutations

find_nonsynon_mutations: find_nonsynon_mutations.hs dirty_hacks.cabal
	cabal install . --bindir=.

test: install
	./find_nonsynon_mutations ./inputs/2.water ./inputs/2.aln

.phony: install test

