

install: find_nonsynon_mutations

find_nonsynon_mutations: find_nonsynon_mutations.hs dirty_hacks.cabal
	cabal install . --bindir=. --disable-documentation

test: install
	./find_nonsynon_mutations ./inputs/2.water ./inputs/2.aln

test2: install
	./find_nonsynon_mutations ./inputs/995.water ./inputs/995.aln

test3: install
	./find_nonsynon_mutations ./inputs/825.water ./inputs/825.aln
#	./find_nonsynon_mutations combined_alignment_files/825.water combined_alignment_files/825.aln

clean:
	rm -rf ./dist/ ./find_nonsynon_mutations

.phony: install test

