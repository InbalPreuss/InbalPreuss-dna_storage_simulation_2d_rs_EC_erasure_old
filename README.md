# DnaStorage

## Usage

1. set values in config.
2. run:
```console
python main.py
``` 

## Test
```console
python -m pytest
```

## Profiling
in test_dna.py make sure that the function in
```console
if __name__ == '__main__':
```
is:
```console
code_profiling()
```
Then, run:
```console
python -m cProfile -s time test_dna.py > temp_file && head -n 10 temp_file
```


## TODO:
2. Make RS for each oligo
6. After synthesising the DNA we should 
    1) Shuffle the data and then use the sort algorithm to sort the oligo according the barcode
    2) Take M(= density of the sequencing) samples each time for the reading   
7. Code that generates different data sizes
10. Make graphs for different number of oligos copies: 20, 100, 1000, 10000. Each graph should use all the M samples
10. Make graphs for different M: 20, 50, 100. Each graph should use 10000 copies for one oligo.    
11. Make tests for RS


## Inbal:

4. Run KB Data and measure the time
5. Run MB Data and measure the time

## Future:

1. Add Multi threading
3. Make RS for all oligos
12. Make RS for 16c3 
13. Make RS for 16c7

## Done:

1. Make sure when adding and removing a letter, of reading a 3 letter sequence that does not exist, we would remove that specific oligo
8. The error graph should be logarithmic graph
9. Remove the exactly 5 bins to distribute the oligos in the synthesis phase - in real life it could be that not all oligos are chosen to In big numbers it probably wouldn't happen). We can cahnge the number of generated oligos to be 100 and not 20.