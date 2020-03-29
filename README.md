# DnaStorage

## TODO:
1. Add Multi threading
2. Make RS for each oligo
3. Make RS for all oligos
4. Run KB Data and measure the time
4. Run MB Data and measure the time
5. After synthesising the DNA we should 
    1) Shuffle the data and then use the sort algorithm to sort the oligo according the barcode
    2) Take M(= density of the sequencing) samples each time for the reading   
6. Code that generates different data sizes
7. The error graph should be logarithmic graph
8. Remove the exactly 5 bins to distribute the oligos in the synthesis phase - in real life it could be that not all oligos are chosen to In big numbers it probably wouldn't happen). We can cahnge the number of generated oligos to be 100 and not 20.
9. Make graphs for different number of oligos copies: 20, 100, 1000, 10000. Each graph should use all the M samples
10. Make graphs for different M: 20, 50, 100. Each graph should use 10000 copies for one oligo.    
11. Make tests for RS
12. Make RS for 16c3 
13. Make RS for 16c7
14. Make sure when adding and removing a letter, of reading a 3 letter sequence that does not exist, we would remove that specific oligo