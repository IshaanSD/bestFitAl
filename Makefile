

.PHONY: clean cmake best_fit_align suffix_array

cmake:
		mkdir -p build
		cmake -S . -B build

best_fit_align:
		cmake --build build --target BestFitExecutable

suffix_array:
		cmake --build build --target SuffixArrayExecutable

clean:
	-rm -r build