import GeneSequencing


class Matrix:
    def __init__(self, seq1: str, seq2: str, isBanded: bool = False):
        self.seq1 = seq1
        self.seq2 = seq2
        self.matrix = dict()
        self.isBanded = isBanded
        self.finalResultPosition: tuple = self.__calculateFinalResultPosition()
        self.__fillMatrix()

    def __calculateFinalResultPosition(self) -> tuple:
        if not self.isBanded:
            return len(self.seq2), len(self.seq1)

        len2 = len(self.seq2)
        len1 = len(self.seq1)

        if abs(len2 - len1) <= GeneSequencing.MAXINDELS:
            return len2, len1

        if len2 < len1:
            return len2, len2 + GeneSequencing.MAXINDELS
        else:
            return len1 + GeneSequencing.MAXINDELS, len1

    def __fillMatrix(self):
        self.__initMatrix()

        for row in range(1, len(self.seq2) + 1):
            leftBound, rightBound = self.__findColumnBounds(row)
            if rightBound - leftBound < 1:  # when there are way more rows than columns
                break

            for col in range(leftBound, rightBound):
                leftPos = (row, col - 1)
                upPos = (row - 1, col)
                diagPos = (row - 1, col - 1)

                finalVal = float("inf")
                prev = (-1, -1)

                if self.matrix.__contains__(leftPos):
                    leftVal = self.matrix[leftPos][0]
                    if leftVal < finalVal:
                        finalVal = self.matrix[leftPos][0] + GeneSequencing.INDEL
                        prev = leftPos
                if self.matrix.__contains__(upPos):
                    upVal = self.matrix[upPos][0] + GeneSequencing.INDEL
                    if upVal < finalVal:
                        finalVal = upVal
                        prev = upPos
                if self.matrix.__contains__(diagPos):
                    diagVal = self.matrix[diagPos][0] + self.__diff(self.seq2[row - 1], self.seq1[col - 1])
                    if diagVal < finalVal:
                        finalVal = diagVal
                        prev = diagPos

                self.__add((row, col), finalVal, prev)

    def __findColumnBounds(self, row: int):
        if not self.isBanded:
            return 1, len(self.seq1) + 1

        leftBound = 1
        temp = row - GeneSequencing.MAXINDELS
        if temp >= 1:
            leftBound = temp
        rightBound = row + GeneSequencing.MAXINDELS + 1
        if rightBound > len(self.seq1) + 1:
            rightBound = len(self.seq1) + 1

        return leftBound, rightBound

    def getEditDistance(self):
        return self.matrix[self.finalResultPosition][0]

    def getEditedSequences(self):
        res1 = ""
        res2 = ""

        currPos = self.finalResultPosition

        while currPos != (0, 0):
            prevPos = self.__getPrev(currPos)
            direction = self.__getDirection(currPos, prevPos)

            if direction == "d":
                res1 += self.seq1[currPos[1] - 1]
                res2 += self.seq2[currPos[0] - 1]

            elif direction == "u":
                res1 += "-"
                res2 += self.seq2[currPos[0] - 1]

            else:
                res1 += self.seq1[currPos[1] - 1]
                res2 += "-"

            currPos = prevPos

        return res1[::-1], res2[::-1]

    def __getPrev(self, pos: tuple) -> tuple:
        return self.matrix[pos][1]

    def __getDirection(self, pos: tuple, prev: tuple) -> str:
        if prev[0] == pos[0] - 1 and prev[1] == pos[1] - 1:
            return "d"  # diagonal

        if prev[0] == pos[0] - 1:
            return "u"  # up

        if prev[1] == pos[1] - 1:
            return "l"  # left

        raise Exception("invalid prev pos")

    def __diff(self, char1: str, char2: str) -> int:
        if char1 == char2:
            return GeneSequencing.MATCH

        return GeneSequencing.SUB

    def __add(self, pos: tuple, val: int, prev: tuple):
        self.matrix[pos] = (val, prev)

    def __initMatrix(self):
        self.__add((0, 0), 0, (-1, -1))

        maxDistance1 = len(self.seq1) + 1
        if self.isBanded:
            maxDistance1 = GeneSequencing.MAXINDELS + 1

        for col in range(1, maxDistance1):
            self.__add((0, col), GeneSequencing.INDEL * col, (0, col - 1))

        maxDistance2 = len(self.seq2) + 1
        if self.isBanded:
            maxDistance2 = GeneSequencing.MAXINDELS + 1

        for row in range(1, maxDistance2):
            self.__add((row, 0), GeneSequencing.INDEL * row, (row - 1, 0))

    def debug_show_matrix(self):
        for row in range(0, len(self.seq2) + 1):
            line = ""
            for col in range(0, len(self.seq1) + 1):
                if self.matrix.__contains__((row, col)):
                    line += str(self.matrix[(row, col)][0]).ljust(7, " ")
                else:
                    line += "_".ljust(7, " ")

            print(line)


if __name__ == '__main__':
    # d = dict()
    # d[(1, 0)] = 1
    # a: bool = d.__contains__((1, 1));
    m = Matrix("acgtacgtacgtacatcgatcgatcgatcagatatatagcgcatagcgatcgtacgtacgcgatatcgtaagatacttctcgagat",
               "satatatacgcgcgtgtgtcacgagctaggagatcgcgatgagcgctcgcgagacacagtatatatatcgagagcgcgatagagta", True)
    m.debug_show_matrix()
    a, b = m.getEditedSequences()
    ed = m.getEditDistance()

    print(a)
    print(b)
    print(ed)
