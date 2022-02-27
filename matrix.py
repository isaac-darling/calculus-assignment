#Isaac Darling, May 2021
#implement several matrix operations and some of their use cases

class Matrix:
    #init Matrix object with row count, column count, and matrix as attributes
    def __init__(self, r, c, vals = None):
        self.r = r
        self.c = c
        self.mat = [[0]*c for _ in range(r)]
        if vals is not None:
            if isinstance(vals[0], list):
                self.mat = list(vals)
            else:
                n = len(vals)
                count = 0
                for i in range(r):
                    for j in range(c):
                        if count >= n:
                            break
                        self.mat[i][j] = vals[count]
                        count+=1
                    else:
                        continue
                    break

    #multiuse function for determinig validity of operations (ensuring their is a matrix with valid dimensions)
    def match(self, other, rev = False):
        if not isinstance(other, Matrix):
            return False

        if rev and self.c != other.r:
            return False

        if not rev and (self.r != other.r or self.c != other.c):
            return False
        return True

    #returns specified minor of a matrix and row/column to ignore
    def minor(self, r, c):
        return Matrix(self.r-1, self.c-1, [row[:c]+row[c+1:] for row in (self.mat[:r]+self.mat[r+1:])])

    #returns Matrix_ji given Matrix_ij
    def transpose(self):
        return Matrix(self.c, self.r, [self.mat[i][j] for j in range(self.c) for i in range(self.r)])

    #computes determinant of matrix (either the input or self if there is no input)
    def determinant(self, matrix = None):
        matrix = matrix or self

        if not isinstance():
            raise TypeError("the determinant is only computable for matrices")

        if matrix.r != matrix.c:
            raise ValueError("the matrix is not square, and therefore has no determinant")

        if matrix.r == 2:
            return (matrix.mat[0][0]*matrix.mat[1][1]) - (matrix.mat[0][1]*matrix.mat[1][0])

        determinant = 0
        for r in range(matrix.r):
            determinant+=((-1)**r)*matrix.mat[0][r]*matrix.determinant(matrix.minor(0, r))
        return determinant

    #returns the multiplicitive inverse of self
    def inverse(self):
        determinant = self.determinant()
        if determinant == 0:
            raise ValueError("the matrix has determinant 0, and therefore no inverse")
        inv_determinant = 1/determinant

        if self.r == 2:
            return inv_determinant*Matrix(2, 2, [self.mat[1][1], -1*self.mat[0][1], -1*self.mat[1][0], self.mat[0][0]])

        return inv_determinant*Matrix(self.r, self.c, [((-1)**(i+j))*self.minor(i, j).determinant() for i in range(self.r) for j in range(self.c)]).transpose()

    #solves linear system of equations using the inverse of the coefficient matrix
    @staticmethod
    def solve(coefficients, variables, constants):
        if not (isinstance(coefficients, Matrix) and isinstance(variables, Matrix) and isinstance(constants, Matrix)):
            raise TypeError("this function only accepts three valid matrices as input")

        prod = coefficients.inverse()*constants
        return dict(zip([row[0] for row in variables.mat], [int(row[0]) if int(row[0]) == row[0] else row[0] for row in prod.mat]))

    #solves linear system of equations using Cramer's rule
    @staticmethod
    def cramer(coefficients, variables, constants):
        if not (isinstance(coefficients, Matrix) and isinstance(variables, Matrix) and isinstance(constants, Matrix)):
            raise TypeError("Cramer's rule must be applied to three valid matrices")

        determinant = coefficients.determinant()
        if determinant == 0:
            raise ValueError("the matrix has determinant 0, and therefore Cramer's rule does not apply")
        inv_determinant = 1/determinant

        determinants = []
        n = len(variables)

        for c in range(n):
            determinants.append(inv_determinant*Matrix(n, n, [constants.mat[i][0] if j == c else coefficients.mat[i][j] for i in range(n) for j in range(n)]).determinant())

        return dict(zip([row[0] for row in variables.mat], [int(x) if int(x) == x else x for x in determinants]))

    #returns the sum of two valid matrices
    def __add__(self, other):
        if not self.match(other):
            raise ArithmeticError("matrices can only be added if their dimensions match")
        return Matrix(self.r, self.c, [self.mat[i][j]+other.mat[i][j] for i in range(self.r) for j in range(self.c)])

    #returns the difference of two valid matrices
    def __sub__(self, other):
        if not self.match(other):
            raise ArithmeticError("matrices can only be subtracted if their dimensions match")
        return Matrix(self.r, self.c, [self.mat[i][j]-other.mat[i][j] for i in range(self.r) for j in range(self.c)])

    #returns the product of two matrices or a matrix and a scalar
    def __mul__(self, other):
        if self.match(other, True):
            vals = []
            for i in range(self.r):
                for j in range(other.c):
                    val = 0
                    for k in range(self.c):
                        val+=self.mat[i][k]*other.mat[k][j]
                    vals.append(val)
            return Matrix(self.r, other.c, vals)

        elif isinstance(other, Matrix):
            raise ValueError("invalid dimensions for matrix multiplication")

        elif isinstance(other, (int, float)):
            return Matrix(self.r, self.c, [self.mat[i][j]*other for i in range(self.r) for j in range(self.c)])

        else:
            raise TypeError("only defined for matrices and real numbers")

    #calls __mul__ but allows for repositioning of the '*' operator
    def __rmul__(self, other):
        return self.__mul__(other)

    #returns whether two matrices are equivalent
    def __eq__(self, other):
        if not isinstance(other, Matrix):
            return False
        return self.mat == other.mat

    #len() returns the total number of elements
    def __len__(self):
        return self.r * self.c

    #outlines the string representation of a Matrix object
    def __repr__(self):
        return str(self.mat).replace("],", "],\n")