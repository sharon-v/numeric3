""" assignment 3
    team members: Hadar Amsalem, Sharon Vazana
    git link: https://github.com/sharon-v/numeric3.git"""







def dominantDiagonal(a):
    sum1 = 0
    for i in range(len(a)):
        for j in range(len(a[i])):
            if i != j:
                sum1 += a[i][j]
        if abs(sum1) > a[i][i]:
            return False
    return True


def guessJ(G, H, b):
    epsilon = 0.00001
    iteration = 0
    prevX = makeMatrics(len(b), len(b[0]))
    currentX = unitMatrics(makeMatrics(len(b), len(b[0])))
    while abs(currentX[0]-prevX[0]) > epsilon:
        for i in range(len(prevX)):
            prevX[i] = plusMatrix(multMatrics(G, prevX[i]), H)
            print(prevX[i], end="   ")
        iteration += 1
        currentX = prevX
        print("")
    print("Total iterations: " + iteration)



def multScalar(a, s):
    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] *= s


def plusMatrix(a, b):
    plusM = unitMatrics(makeMatrics(len(a), len(a[0])))
    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] += b[i][j]
    return plusM


def minusMatrix(a, b):
    minusM = unitMatrics(makeMatrics(len(a), len(a[0])))
    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] -= b[i][j]
    return minusM