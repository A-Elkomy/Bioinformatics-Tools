import numpy as np
#Ahmed Tarek Elkomy
#Mostafa mohamed mahmoud
#Marwan Akram Said
#Moahmed Afifi Abd el hafez
#MMohamed Ahmed Moahmed Moastafa
if __name__ == '__main__':
    dtn = [[0, 20, 60, 100, 90],
           [20, 0, 50, 90, 80],
           [60, 50, 0, 40, 50],
           [100, 90, 40, 0, 30],
           [90, 80, 50, 30, 0]]
    dtn = np.array(dtn)
    print(dtn)
    print()
    print()
    Sequence_names = ["A", "B", "C", "D", "E"]
    seq_num = len(Sequence_names)
    while (seq_num > 2):
        minNumber = 888
        for i in range(len(dtn)):
            for j in range(len(dtn)):
                if dtn[i][j] == 0:
                    continue
                else:
                    if minNumber > dtn[i][j]:
                        minNumber = dtn[i][j]

        print(f"Min = {minNumber}")
        print(dtn)
        print('-' * 50)
        indexMinNumber = list()
        index_x, index_y = np.where(dtn == minNumber)
        for i in range(len(index_x)):
            indexMinNumber.append([index_x[i], index_y[i]])
        row = indexMinNumber[0][0]
        col = indexMinNumber[0][1]
        seq1 = -1
        seq2 = -1
        if row > col:
            seq1 = col
            seq2 = row
        else:
            seq1 = row
            seq2 = col
        for i in range(0, seq1):
            dtn[i][seq1] = (dtn[i][seq1] + dtn[i][seq2]) / 2
            dtn[seq1][i] = dtn[i][seq1]

        for i in range(seq1 + 1, seq_num):
            if i < seq2:
                dtn[seq1][i] = (dtn[seq1][i] + dtn[i][seq2]) / 2
                dtn[i][seq1] = dtn[seq1][i]
            else:
                dtn[seq1][i] = (dtn[seq1][i] + dtn[seq2][i]) / 2
                dtn[i][seq1] = dtn[seq1][i]
        dtn = np.delete(dtn, seq2, axis=0)
        dtn = np.delete(dtn, seq2, axis=1)
        Sequence_names[seq1] = ("(" + Sequence_names[seq1] + "," + Sequence_names[seq2] + ")")
        del Sequence_names[seq2]
        seq_num = seq_num - 1
        if (seq_num == 2):
            print(f"Min = {dtn[0][1]}")
            print(dtn)
            print('-' * 50)
            Sequence_names[0] = ("(" + Sequence_names[0] + "," + Sequence_names[1] + ")")
            del Sequence_names[1]
    print(Sequence_names)
