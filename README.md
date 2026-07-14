# LLL Subset-Sum Experiments

This repository contains a Fortran implementation of experimental lattice-based attacks on the **subset-sum problem** using the **Lenstra–Lenstra–Lovász (LLL) lattice basis reduction algorithm**.

The program generates random subset-sum instances, embeds each instance into a lattice, applies LLL reduction, and searches the reduced basis for a vector encoding a valid subset-sum solution.

## The Subset-Sum Problem

Given positive integers

$$
s_1,s_2,\ldots,s_n
$$

and a target integer $t$, the subset-sum problem asks whether there exists a binary vector

$$
x=(x_1,x_2,\ldots,x_n), \qquad x_i\in{0,1},
$$

such that

$$
\sum_{i=1}^{n}s_i x_i=t.
$$

In each experiment, the program first generates a random hidden binary vector $y$ and constructs the target

$$
t=\sum_{i=1}^{n}s_i y_i.
$$

Thus, every generated instance is guaranteed to have at least one solution.

## Lattice Approach

The program constructs a lattice basis associated with the subset-sum instance and applies LLL basis reduction. The purpose of LLL reduction is to find relatively short vectors in the lattice.

For the embedding used in this program, a vector corresponding to a subset-sum solution has the form

$$
\left(x_1-\frac12,\ldots,x_n-\frac12,0\right),
$$

up to sign, where each $x_i$ is either 0 or 1. Consequently, the first $n$ coordinates are all either $+1/2$ or $-1/2$, while the final coordinate is zero.

After LLL reduction, the program examines the reduced basis columns for vectors having this form. Every candidate is verified directly by checking whether

$$
\sum_{i=1}^{n}s_i x_i=t.
$$

## What the Program Does

For each randomly generated problem, the program:

1. Generates $n$ positive integers $s_1,\ldots,s_n$.

2. Generates a random hidden binary subset.

3. Computes the corresponding target sum.

4. Calculates the approximate subset-sum density

   $$
   d=\frac{n}{\log_2(\max_i s_i)}.
   $$

5. Constructs a lattice basis for the subset-sum instance.

6. Applies LLL lattice basis reduction.

7. Searches the reduced basis columns for an encoded subset-sum solution.

8. Verifies any recovered solution directly.

9. Reports whether the recovered subset is identical to the originally generated hidden subset.

10. Displays an overall success rate after all experiments are complete.

## Implementation

The source code is contained in:

```text
experiments.f90
```

The program contains:

* A Gram–Schmidt orthogonalization routine.
* An LLL size-reduction routine.
* An LLL basis-reduction routine.
* A routine for identifying and verifying subset-sum solution vectors.
* User-friendly reporting of the hidden subset, target, reduced basis column norms, recovered solution, and overall experimental success rate.

Basis vectors are stored as **columns** of the lattice basis matrix.

The implementation uses double-precision floating-point arithmetic for lattice computations and 64-bit integers for subset-sum values and indices.

For simplicity and reliability in small experiments, the Gram–Schmidt data are recomputed after each modification of the lattice basis. This is slower than a highly optimized LLL implementation, but it makes the code easier to understand and experiment with.

## Default Parameters

The current version uses:

```fortran
integer(int64), parameter :: n = 20
integer(int64), parameter :: number_of_problems = 10
real(dp), parameter :: delta = 0.85_dp
```

Thus, by default, the program generates ten random subset-sum problems, each containing twenty integers.

These values can be changed directly in `experiments.f90`.

## Compiling

A Fortran compiler supporting the `iso_fortran_env` intrinsic module is required.

Using `gfortran`:

```bash
gfortran -Wall -Wextra -O2 experiments.f90 -o experiments
```

### Windows

Run with:

```text
experiments.exe
```

### Linux or macOS

Run with:

```bash
./experiments
```

## Output

For each experiment, the program displays:

* The generated subset-sum integers.
* The hidden binary subset.
* The values selected by the hidden subset.
* The target sum.
* The approximate density of the instance.
* The norms of the LLL-reduced basis columns.
* Whether a solution was recovered.
* The recovered binary subset, if found.
* Direct verification that the recovered subset produces the target.
* Whether the recovered subset is identical to the original hidden subset.

At the end, the program reports the total number of problems tested, the number of solutions recovered, and the overall success rate.

## Experimental Nature

This code is intended for experimentation and education rather than as an optimized cryptographic lattice library.

In particular:

* The LLL implementation uses floating-point arithmetic.
* Gram–Schmidt data are recomputed frequently for simplicity.
* The solution search examines the reduced basis columns directly.
* LLL does not guarantee that a desired short vector will appear as one of the reduced basis vectors.
* Failure to recover a solution does not mean that the subset-sum instance has no solution; every generated instance has a known solution by construction.

The program is therefore useful for studying experimentally when lattice reduction exposes the hidden structure of subset-sum instances.

## Historical Background

The use of lattice basis reduction to attack low-density subset-sum problems is associated with work by Andrew M. Odlyzko and Jeffrey C. Lagarias. Their work demonstrated that sufficiently low-density subset-sum problems can often be attacked by transforming them into problems involving unusually short vectors in lattices.

The LLL algorithm, introduced by Arjen Lenstra, Hendrik Lenstra, and László Lovász, provides a polynomial-time method for finding reduced lattice bases containing relatively short vectors.

This repository is an experimental implementation inspired by these ideas.

## Author

**Craig Alan Feinstein**

## License and Provenance

This repository contains an experimental Fortran implementation of LLL lattice reduction and its application to subset-sum problems.

The program was developed from an earlier Fortran LLL implementation found online and subsequently rewritten, corrected, expanded, and documented. The original online source has not yet been identified.

Until the provenance and licensing terms of the earlier implementation are established, no open-source license is granted for this repository.
