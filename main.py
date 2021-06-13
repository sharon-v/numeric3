""" assignment 3
    team members: Hadar Amsalem, Sharon Vazana
    git link: https://github.com/sharon-v/numeric3.git"""

# Gauss–Seidel
# Jacobi

# ########## LDU method ############
from _json import make_encoder


def makeMatrics(row, col=1):
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
    temp = c[r1]
    c[r1] = c[r2]
    c[r2] = temp
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
        print("iter. " + str(iteration) + " =", end="\t")
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


# ########### end LDU method #############

##########################################

# ########## iterative method ############
def jacobiIter(a, b):
    """
    :param a: coefficient matrics
    :param b: solution vector
    :return: prints Jacobi method iterations
    """
    epsilon = 0.00001
    iteration = 0
    flag = True
    prevX = makeMatrics(len(a))  # start as zero vector
    currentX = makeMatrics(len(a))
    matA, vectorB = isolateVariables(a, b)
    print("\n*** iter Jacobi *** ")
    print("-----------------")
    while abs(currentX[0][0] - prevX[0][0]) > epsilon or flag is True:
        flag = False
        prevX = currentX
        currentX = makeMatrics(len(currentX))
        if iteration >= 100:
            print("The system can't converge :(")
            break
        print("iter. " + str(iteration) + " =", end="\t")
        print(prevX)
        for i in range(len(a)):
            j = 0
            currentX[i][0] = vectorB[i][0]
            while j < len(a[0]):
                if j is not i:
                    currentX[i][0] += matA[i][j] * prevX[j][0]
                j += 1
        iteration += 1
    print("-----------------")
    print("Total number of iterations: " + str(iteration))
    print("-----------------")


def gaussSeidelIter(a, b):
    """
    :param a: coefficient matrics
    :param b: solution vector
    :return: prints Gauss Seidel method iterations
    """
    epsilon = 0.00001
    iteration = 0
    flag = True
    prevX = makeMatrics(len(a))  # start as zero vector
    currentX = makeMatrics(len(a))
    matA, vectorB = isolateVariables(a, b)
    print("\n*** iter Gauss-Seidel *** ")
    print("-----------------")
    while abs(currentX[0][0] - prevX[0][0]) > epsilon or flag is True:
        flag = False
        prevX[0][0] = currentX[0][0]
        if iteration >= 100:
            print("The system can't converge :(")
            break
        print("iter. " + str(iteration) + " =", end="\t")
        print(currentX)
        for i in range(len(a)):
            j = 0
            currentX[i][0] = vectorB[i][0]
            while j < len(a[0]):
                if j is not i:
                    currentX[i][0] += matA[i][j] * currentX[j][0]
                j += 1
        iteration += 1
    print("-----------------")
    print("Total number of iterations: " + str(iteration))
    print("-----------------")


def isolateVariables(a, b):
    """
    :param a: coefficient matrics
    :param b: solution vector
    :return: new matrics & vector after isolation of pivot in each row
    """
    vectorB = makeMatrics(len(a))
    matA = makeMatrics(len(a), len(a[0]))
    for i in range(len(a)):
        vectorB[i][0] = b[i][0] / a[i][i]
        j = 0
        while j < len(a[0]):
            if j is i:
                matA[i][i] = 1
            else:
                matA[i][j] -= a[i][j] / a[i][i]
            j += 1
    return matA, vectorB


# ######### end iterative method ##########
###########################################
# ########## dominant diagonal ############
def copyMat(A):
    """
    :param A: a matrix
    :return: a copy of the matrix
    """
    B = makeMatrics(len(A), len(A[0]))  # create a zero matrix of the same size as A
    for i in range(len(A)):  # copy A
        for j in range(len(A[0])):
            B[i][j] = A[i][j]
    return B


def createDominantDiagonal(A, b=None):
    """
    :param A: a coefficients matrix
    :param b: the column vector of constant terms.
    :return: matrix A with dominant diagonal
    """
    max = 0  # the max value in the current row or column in the matrix
    maxIndex = 0  # the index of the max value
    for i in range((len(A))):  # calc the max value for each member on the main diagonal
        sum = 0  # the sum of the members in the current row in A
        for j in range(len(A)):  # go over the current row
            sum += abs(A[i][j])  # add the value of each member in the row
            if abs(A[i][j]) > max:  # search for the max value in the current row
                max = abs(A[i][j])
                maxIndex = j
        if (sum - max) <= max:  # if the max value in the row meets the condition of a dominant diagonal
            A = manualSwapCol(A, maxIndex, i)  # swap between the columns of the ...
            # ... current value on the main diagonal and the max value in that row
        else:  # look for the max value in the current column
            max = 0
            maxIndex = 0
            for j in range(len(A)):  # go over the current column
                # sum += abs(A[j][i])
                if abs(A[j][i]) > max:  # search for the max value in the current column
                    max = abs(A[j][i])
                    maxIndex = j
            if rowSum(A[j]) - max <= max:  # if the max value in the row meets the condition of a dominant diagonal
                A, b = manualSwapRow(A, b, i, maxIndex)  # swap between the rows of the current value on...
                # ...the main diagonal and the max value in that column
            else:
                print("ERROR - no dominant diagonal")  # A can't be changed into dominant diagonal matrix
                return None, None
    return A, b


def manualSwapRow(a, b, r1, r2):
    """
    manaul rows exchange (without e)
    :param a: The coefficient matrix
    :param b:  The column vector of constant terms
    :param r1: the first row to swap
    :param r2: the second row to swap
    :return: the matrix after the swap, The column vector of constant terms after swap
    """
    if r2 < len(a) and r1 < len(a):
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp
        if b is not None:  # if the result vector is not none swap him too
            temp = b[r1]
            b[r1] = b[r2]
            b[r2] = temp
    return a, b


def manualSwapCol(a, c1, c2):
    """
    :param a: The coefficient matrix
    :param c1: the first column to swap
    :param c2: the second column to swap
    :return: the matrix after the swap
    """
    if c2 < len(a) and c1 < len(a):
        for i in range(len(a)):
            temp = a[i][c1]
            a[i][c1] = a[i][c2]
            a[i][c2] = temp
    return a


def rowSum(line):
    """
    :param line: A list od numbers - line for the matrix
    :return: the sum of all the numbers in abs  in the list
    """
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


# dominant diagonal end


def driver():
    """
    main function for GH method
    :return: prints results
    """
    a = [[4, 2, 0],
         [2, 10, 4],
         [0, 4, 5]]

    b = [[2],
         [6],
         [5]]

    if dominantDiagonal(a) is False:
        # check
        copyA = copyMat(a)
        copyB = copyMat(b)
        copyA, copyB = createDominantDiagonal(copyA, copyB)
        if (copyA is not None) and (copyB is not None):
            a = copyA
            b = copyB
        else:
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


def newDrive():
    """
    main function for iterative isolation of variables method
    :return: prints results
    """
    a = [[4, 2, 0],
         [2, 4, 10],
         [0, 4, 5]]

    b = [[2],
         [6],
         [5]]

    if dominantDiagonal(a) is False:
        # check
        copyA = copyMat(a)
        copyB = copyMat(b)
        copyA, copyB = createDominantDiagonal(copyA, copyB)
        if (copyA is not None) and (copyB is not None):
            a = copyA
            b = copyB
        else:
            print("No dominant diagonal")
    print("enter 0 for jacobi, 1 for gauss-seidel, other for both: ", end="")
    temp = input()
    if temp is '0':
        jacobiIter(a, b)
    elif temp is '1':
        gaussSeidelIter(a, b)
    else:
        jacobiIter(a, b)
        gaussSeidelIter(a, b)


newDrive()
