# math_tools
Pretty much what it says: Math Tools

## To Do
- Modularize no-extension files
- Documentation
- Integrate Vectors and PolyVectors into the Tensor namespace
  
## Differential
Includes a function that numerically approximates the derivative of a function at a discrete point x.

Also includes functions that numerically approximate the roots of a given function with various methods.

### Functions in `differential`:
  ```
  diff(*function_pointer, x_point [default 0], dx [default 10^-5])
  ```
  Finds the derivative of the function where `*function_pointer` is a pointer to a function that takes in a single argument of type `double` and returns a `double`;
  `x_point` is the x-value you want to evaluate the derivative at; and `dx` is the optional step parameter - too big and your
  derivative is inaccurate, too small and it goes haywire to + or - infinity.

  ```
  bisectRoot(*function_pointer, interval{1, 2}, steps)
  ```
  Approximates the roots via the bisection method. Same requirements for `*function_pointer` as above. `interval` is an `std::array`
  of the interval you want the function to work on *(I just insert `{1, 2}` lists and let the STL initialize the array with that)*.
  `steps` is the maximum number of recursions this function can perform before returning the approximated result.

  Works best when the function results to opposite signs when evaluated at `interval[0]` and `interval[1]`.

  ```
  newtonRoot(*function_pointer, initial_guess, steps)
  ```
  Approximates roots via Newton's method. Same function requirements. Steps are the used in the same way as above. `initial_guess` is
  the initial value you want the function to work on.

  Works best when there are no non-zero local extrema close to the function evaluated at `initial_guess`.

  ```
  secantRoot(*function_pointer, interval{1, 2}, steps)
  ```
  Roots via the method of secants. Same description and use as `bisectRoot`.

## Integral
Methods to numerically approximate definite integrals and solve for differential equations. All functions are in the namespace `integral`.
Can evaluate polynomials too if included after `#polynomial`.

  ### Functions in `integral`:
  ```
  integral::lRect(*function_pointer, start, stop, terms [default 4096])
  ```
  Approximates via the left rectangle Riemann sum. Same requirements for the function. `start` is the lower
  bound, and `stop` is the upper bound. `terms` is the number of terms in the summation to approximate the
  integral: higher = better approximation, but uses more resources.
  ```
  integral::mRect(*function_pointer, start, stop, terms [default 4096]);
  integral::rRect(*function_pointer, start, stop, terms [default 4096]);
  integral::trapInt(*function_pointer, start, stop, terms [default 4096]);
  integral::simpsonInt(*function_pointer, start, stop, terms [default 4096]);
  ```
  The above functions have basically the same description, only differing in methods used:
  - `mRect(...)` uses the midpoint rectangle.
  - `rRect(...)` uses the righ rectangle.
  - `trapInt(...)` uses the trapezoid method: can be more or less accurate depending on the function. Uses more resources.
  - `simpsonInt(...)` is a step further from the trapezoid method. Approximates with parabolic shapes instead. Uses more resources than `trapInt(...)`.

  ### Differential Equations
  I included them here, under the namespace `diffeqn`.
  - `diffeqn::euler(pos, vel, time, *Dp, *Dv)` - uses the modified Euler integration method.
  - `diffeqn::sEuler(pos, vel, time, *Dp, *Dv)` - uses the simple Euler integration method. Less resources, less accuracy.
  - `diffeqn::rkm(pos, vel, time, *Dp, *Dv)` - uses the Runge-Kutta (RK4) method. More accurate, more resource use.
  ### WARNING: These do not return any value. They modify the original p, v, and t parameters. You will have to manually record the values after each iteration of these functions.
  Parameter descriptions:
  - `pos`, `vel`, `time` - current position, velocity, and time. Must be of type `double`
  - `*Dp` and `*Dv` - pointers to the instantaneous velocity and instantaneous acceleration functions, respectively. Each of these functions
    must take three `double` parameters in the order `(current_position, current_velocity, current_time)`.

  After passing the appropriate arguments, the function *updates* the initial values passed into it using the given velocity and acceleration functions
  with the internally-defined time step `dt`. Again, you must make a record of the values per iteration otherwise the previous values are lost to time.

## Interpolate (Problems Exist)
Interpolate to find polynomial functions for a given 2D dataset. Since there are still problems, I will not write the documentation for now.

## Linear
Helper functions. You might use them, you might not. Include when necessesary.
- `linspace(start, stop, terms)` gives you a vector containing `terms` number of evenly spaced `double` terms
  between `start` and `stop` (inclusive).
- `abs(x)` returns the absolute value of `x`.
- `signof(x)` returns an `std::string` containing the sign of `x` (either `+` or `-`).

## Matrix
Matrix stuff in the namespace `Tensor`. Stores matrix data as a vector of vectors of doubles.
### Class
```
Tensor::Matrix {
  private:
  uint32_t length, width;
  std::vector<std::vector<double>> DATA;

  public:
  // Constructors and Member Functions
};
```

### Constructors
```
Matrix(W, H);
```
Constructs a matrix with width `W` and height `H`. The matrix is filled with 0.0s of the type `double`.

```
Matrix(int N);
```
Constructs an `N`x`N` unit matrix *(i.e., 1.0s on the upper left-lower right diagonal, zeroes everywhere else)*.

```
Matrix(std::vector<std::vector<double>>);
```
Constructs a matrix based on the initializer list. The list should be in the form of:
```
std::vector{   //Rows
  std::vector{ // Row 1 column 1, Row 1 column 2, ... }
  std::vector{ // Row 2 column 1, Row 2 column 2, ... }
  ...
}
```

### Member Functions and Operators
```
Matrix::rows();
Matrix::columns();
Matrix::data();
```
Returns a copy of the data they represent.
- `rows()` - the height of the matrix.
- `columns()` - the width of the matrix.
- `data()` returns the entire `DATA` vector of vectors.

```
Matrix::set(row, column, number);
```
Change the data of the matrix at `[row][column]` to `number`.

```
Matrix::determinant();
```
Returns the determinant of a square matrix; will throw an exception and return 0.0 when used on a non-square matrix.

```
GaussSolution(Matrix, std::vector<double>)
```
Finds the solution to a system of equations that can be expressed as an augmented square matrix using Gaussian elimination.

```
Matrix::eighenvalues();
```
Returns the real eigenvalues of a square matrix. Requires `polynomial` to be included before `matrix`.

```
std::cout << Matrix;   // and iostream/fstream variants
```
Prints the matrix's values on screen (or in a file).

### Primitive Operator Overloads
- `+` - matrix addition in the mathematical context.
- `-` - matrix subtraction in the mathematical context.
- `*` - matrix multiplication. **NOT YET IMPLEMENTED.**

```
Tensor::PolyMatrix {}
```
Essentially the same as `Tensor::Matrix`, but uses `polynomial` instead of `double` as its main data type. Has the same
methods implemented *(except `GaussSolution` and `eigenvalues`)*.

## Polynomial
Under construction
