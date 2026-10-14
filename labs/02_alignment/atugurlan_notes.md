## Datasets used

I used three TP53 records: `NG_017013.2`, `PZ385140.1`, and `PZ385139.1`.

I compared two pairs: indices 0 and 1, where one sequence is significantly longer than the other, and indices 1 and 2, where the sequences have the same length and are identical.

## When is global alignment preferred over local alignment?

Global alignment is preferred when we want to see how two complete sequences differ. It compares them from beginning to end, leaving no region out of the comparison.

Local alignment is preferred when we are interested only in a similar region shared by two sequences. It finds the highest-scoring pair of regions and leaves the rest out of the comparison.
