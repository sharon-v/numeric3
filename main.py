""" assignment 3
    team members: Hadar Amsalem, Sharon Vazana
    git link: https://github.com/sharon-v/numeric3.git"""


# Gauss–Seidel
# Jacobi


def makeMatrics(row, col):
    """
    :param row: get rows of matrix
    :param col: get columns of matrix
    :return: return zero matrix
    """
    c = []
    for i in range(row):
        c += [[0] * col]
    return c


def multMatrics(a, b):
    """
    :param a: get matrix
    :param b: get matrix
    :return: new matrix of mul between them
    """
    if len(a[0]) is len(b):
        c = makeMatrics(len(a), len(b[0]))
        for row in range(len(a)):
            for col in range(len(b[0])):
                for x in range(len(a)):
                    c[row][col] += (a[row][x] * b[x][col])
        return c
    return None


def det(a):
    """
    :param a: matrics
    :return: determinant of matrics a
    """
    if len(a) and len(a[0]) is 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]
    sum1 = 0
    for j in range(len(a[0])):
        if j % 2 == 0:
            sign = 1
        else:
            sign = -1
        sum1 += sign * a[0][j] * det(minor(a, 0, j))
    return sum1


def minor(b, row, col):
    """
    :param b: matrics
    :param row: row to remove
    :param col: column to remove
    :return: matrics b without row and col
    """
    if row >= len(b) and col >= len(b):
        return b
    c = makeMatrics(len(b) - 1, len(b) - 1)
    x = 0
    y = 0
    for i in range(len(b)):
        for j in range(len(b[0])):
            if i is not row and j is not col:
                c[x][y] = b[i][j]
                if y is len(c[0]) - 1:
                    x += 1
                    y = 0
                else:
                    y += 1
    return c


def elementalMatrics(a, i, j):
    """
    :param a: matrics
    :param i: row index
    :param j: column index
    :return: elemental matrics to make a[i][j] = 0
    """
    c = makeMatrics(len(a), len(a[0]))
    c = unitMatrics(c)
    c[i][j] = -1 * (a[i][j] / a[j][j])
    return c


def unitMatrics(c):
    """
    :param c: get matrix
    :return: make her to unit matrix
    """
    for x in range(len(c[0])):
        c[x][x] = 1  # make c a unit matrix
    return c


def swapRow(a, r1, r2):
    """
    :param a: original matrics
    :param r1: 1st row
    :param r2: row to swap with
    :return: elemental swap matrics
    """
    c = makeMatrics(len(a), len(a[0]))
    c = unitMatrics(c)
    c[r1] = a[r2]
    c[r2] = a[r1]
    return c


def findU(a):
    """
    :param a: matrics
    :return: U and inverse L matrices
    """
    invL = unitMatrics(makeMatrics(len(a), len(a[0])))
    for row in range(len(a)):
        j = row + 1
        while j < len(a):
            if a[row][row] is 0:  # LU
                k = j + 1
                while k < len(a):
                    if a[k][row] is not 0:
                        b = swapRow(a, row, k)
                        a = multMatrics(b, a)
                        invL = multMatrics(b, invL)
                        break
                    k += 1
            b = elementalMatrics(a, j, row)
            a = multMatrics(b, a)
            invL = multMatrics(b, invL)
            j += 1
    return a, invL


def oneOnDiagonal(a, matInverse):
    """
    :param a: matrics
    :param matInverse: matrics composed from all the elemental operations on matrics a
    :return: a with 1 on the main diagonal, matInverse updated with the new elemental operations
    """
    b = makeMatrics(len(a), len(a[0]))
    b = unitMatrics(b)
    for i in range(len(a[0])):
        b = unitMatrics(b)
        b[i][i] = 1 / a[i][i]
        a = multMatrics(b, a)
        matInverse = multMatrics(b, matInverse)
    return a, matInverse


def inverse(a):
    """
    :param a: matrics
    :return: inverse matrics a
    """
    if det(a) is 0:
        return
    a, matInverse = findU(a)
    a, matInverse = oneOnDiagonal(a, matInverse)
    size = len(a[0]) - 1
    while size > 0:
        for i in range(size):
            b = elementalMatrics(a, i, size)
            a = multMatrics(b, a)
            matInverse = multMatrics(b, matInverse)
        size -= 1
    return matInverse


def LDU(a):
    """
    :param a: get a matrix
    :return: L, U, D matrix
    """
    L = makeMatrics(len(a), len(a[0]))
    D = makeMatrics(len(a), len(a[0]))
    U = makeMatrics(len(a), len(a[0]))
    for i in range(len(a)):
        for j in range(len(a[0])):
            if i is j:
                D[i][j] = a[i][j]
            elif i < j:
                U[i][j] = a[i][j]
            else:
                L[i][j] = a[i][j]
    return L, D, U


def GHjacobi(a):
    """
    :param a: get a matrix
    :return: H, G matrix of jacobi
    """
    L, D, U = LDU(a)
    invD = inverse(D)
    G = multScalar(multMatrics(invD, plusMatrix(L, U)), -1)
    return G, invD


def GHgauss(a):
    """
    :param a: get a matrix
    :return: H, G matrix of gauss
    """
    L, D, U = LDU(a)
    invLminusD = inverse(plusMatrix(L, D))
    invLplusD = inverse(plusMatrix(L, D))
    G = multScalar(multMatrics(invLminusD, U), -1)
    return G, invLplusD


def dominantDiagonal(a):
    """
    :param a: get a matrix
    :return: true if the matrix have dominant diagonal, else false
    """
    for i in range(len(a)):
        sum1 = 0
        for j in range(len(a[0])):
            if i is not j:
                sum1 += abs(a[i][j])
        if sum1 > a[i][i]:
            return False
    return True


def guess(G, H, b):
    """
    :param G: get G matrix
    :param H: get H matrix
    :param b: get result vector
    :return: print iteration of calculus
    """
    epsilon = 0.00001
    iteration = 0
    prevX = makeMatrics(len(b), len(b[0]))
    currentX = makeMatrics(len(b), len(b[0]))
    flag = True
    print("-----------------")
    while abs(currentX[0][0] - prevX[0][0]) > epsilon or flag is True:
        flag = False
        print(prevX)
        currentX = prevX
        prevX = plusMatrix(multMatrics(G, prevX), multMatrics(H, b))
        iteration += 1
    print("-----------------")
    print("Total number of iterations: " + str(iteration))
    print("-----------------")


def multScalar(a, s):
    """
    :param a: get a matrix
    :param s: get scalar number
    :return: new matrix of scalar mult a
    """
    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] *= s
    return a


def plusMatrix(a, b):
    """
    :param a: get a matrix
    :param b: get b matrix
    :return: new matrix of plus between a and b
    """
    plusM = makeMatrics(len(a), len(a[0]))
    for i in range(len(a)):
        for j in range(len(a[0])):
            plusM[i][j] = a[i][j] + b[i][j]
    return plusM


def infNorm(a):
    """
    :param a: matrics
    :return: infinity norm of matrics a
    """
    norm = 0
    for i in range(len(a[0])):
        sumRow = 0
        for j in range(len(a)):
            sumRow += abs(a[i][j])
        norm = max(sumRow, norm)
    return norm


def checkConvergence(G):
    """
    :param G: a matrics
    :return: true if convergence condition is met
    """
    if infNorm(G) < 1:
        return True
    return False


def gaussSeidel(a, b):
    """
    :param a: matrix a
    :param b: result vector
    :return: print the iteration of gauss seidel
    """
    G, H = GHgauss(a)
    print("\n*** Gauss-Seidel *** ")
    if checkConvergence(G) is False:
        print("The system can't converge :(")
        return
    guess(G, H, b)


def jacobi(a, b):
    """
    :param a: matrix a
    :param b: result vector
    :return: print the iteration of jacobi
    """
    G, H = GHjacobi(a)
    print("\n*** Jacobi *** ")
    if checkConvergence(G) is False:
        print("The system can't converge :(")
        return
    guess(G, H, b)


def driver():
    """
    main function
    :return: prints results
    """
    a = [[4, 2, 0],
         [2, 10, 4],
         [0, 4, 3]]

    b = [[2],
         [6],
         [5]]

    if dominantDiagonal(a) is False:
        print("No dominant diagonal")
    print("enter 0 for jacobi, 1 for gauss-seidel, other for both: ", end="")
    temp = input()
    if temp is '0':
        jacobi(a, b)
    elif temp is '1':
        gaussSeidel(a, b)
    else:
        jacobi(a, b)
        gaussSeidel(a, b)


driver()
