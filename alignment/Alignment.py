import numpy as np
#Ahmed Tarek Elkomy
#Mostafa mohamed mahmoud
#Marwan Akram Said
#Moahmed Afifi Abd el hafez
#MMohamed Ahmed Moahmed Moastafa

def globalAlignmentForDNA():
    firstSequence = input("Enter Sequence 1: ").upper()
    secondSequence = input("Enter Sequence 2: ").upper()
    lengthFirstSequence = len(firstSequence)
    lengthSecondSequence = len(secondSequence)

    # setup Matrix
    mainMatrix = np.zeros((lengthFirstSequence + 1, lengthSecondSequence + 1)).astype(int)
    match_checker_matrix = np.zeros((lengthFirstSequence, lengthSecondSequence)).astype(int)

    # scoring for Alignment
    matchPenalty = 1
    mismatchPenalty = -2
    gapPenalty = -1
    print(f"Penalties => Match = {matchPenalty} , Mismatch = {mismatchPenalty} , Gap = {gapPenalty}")
    # Fill the match checker matrix according to match or mismatch
    for i in range(lengthFirstSequence):
        for j in range(lengthSecondSequence):
            if firstSequence[i] == secondSequence[j]:
                match_checker_matrix[i][j] = matchPenalty
            else:
                match_checker_matrix[i][j] = mismatchPenalty
    # initialization main Matrix
    for i in range(lengthFirstSequence + 1):
        mainMatrix[i][0] = i * gapPenalty
    for j in range(lengthSecondSequence + 1):
        mainMatrix[0][j] = j * gapPenalty

    print('-' * 50)
    print(f'Setup Matrix and Initialization')
    print(mainMatrix)

    # make Needleman-water algorithm for filling matrix

    for idx in range(1, lengthFirstSequence + 1):
        for jdx in range(1, lengthSecondSequence + 1):
            if firstSequence[idx - 1] == secondSequence[jdx - 1]:
                mainMatrix[idx][jdx] = max(mainMatrix[idx - 1][jdx - 1] + matchPenalty,
                                           mainMatrix[idx - 1][jdx] + gapPenalty, mainMatrix[idx][jdx - 1] + gapPenalty)
            elif firstSequence[idx - 1] != secondSequence[jdx - 1]:
                mainMatrix[idx][jdx] = max(mainMatrix[idx - 1][jdx - 1] + mismatchPenalty,
                                           mainMatrix[idx - 1][jdx] + gapPenalty, mainMatrix[idx][jdx - 1] + gapPenalty)

    print('-' * 50)
    print(f'Filling matrix')
    print(mainMatrix)

    # Trace Back
    iTraceBack = lengthFirstSequence
    jTraceBack = lengthSecondSequence
    firstAlignment = ''
    secondAlignment = ''
    while iTraceBack > 0 and jTraceBack >= 0:
        if (iTraceBack >= 0 and jTraceBack >= 0 and
                mainMatrix[iTraceBack][jTraceBack] == mainMatrix[iTraceBack - 1][jTraceBack - 1] +
                match_checker_matrix[iTraceBack - 1][jTraceBack - 1]):
            firstAlignment = firstSequence[iTraceBack - 1] + firstAlignment
            secondAlignment = secondSequence[jTraceBack - 1] + secondAlignment
            iTraceBack = iTraceBack - 1
            jTraceBack = jTraceBack - 1
        elif (iTraceBack >= 0 and
              mainMatrix[iTraceBack][jTraceBack] == mainMatrix[iTraceBack - 1][jTraceBack] + gapPenalty):
            firstAlignment = firstSequence[iTraceBack - 1] + firstAlignment
            secondAlignment = '-' + secondAlignment
            iTraceBack = iTraceBack - 1
        else:
            firstAlignment = '-' + firstAlignment
            secondAlignment = secondSequence[jTraceBack - 1] + secondAlignment
            jTraceBack = jTraceBack - 1

    print('-' * 50)
    print(f'Global Alignment and Score')
    print(f"Align1: {firstAlignment}")
    print(f"Align2: {secondAlignment}")
    score = mainMatrix[lengthFirstSequence][lengthSecondSequence]
    print(f"Score = {score}")


def localAlignmentForDNA():
    S1 = input("Enter Longest sequence: ")
    S2 = input("Enter Shortest sequence: ")

    lengthFirstSequence = len(S1)
    lengthSecondSequence = len(S2)
    Matrix = np.zeros((lengthSecondSequence + 1, lengthFirstSequence + 1), dtype=np.int64)

    match = 1
    mismatch = -2
    gap = -1
    print(f"Penalties => Match = {match} , Mismatch = {mismatch} , Gap = {gap}")
    print('-' * 40)

    highmax = -999
    maxarr = []
    for i in range(1, lengthSecondSequence + 1):
        for j in range(1, lengthFirstSequence + 1):
            if S2[i - 1] == S1[j - 1]:
                Matrix[i][j] = max(Matrix[i - 1][j - 1] + match, Matrix[i - 1][j] + gap, Matrix[i][j - 1] + gap, 0)
                if Matrix[i][j] >= highmax:
                    highmax = Matrix[i][j]
            else:
                Matrix[i][j] = max(Matrix[i - 1][j - 1] + mismatch, Matrix[i - 1][j] + gap, Matrix[i][j - 1] + gap, 0)
                if Matrix[i][j] >= highmax:
                    highmax = Matrix[i][j]

    print("Filled Matrix")
    print(Matrix)
    print('-' * 40)

    print(f'Max Alignment score = {highmax}')
    print('-' * 40)

    index_x, index_y = np.where(Matrix == highmax)
    for i in range(len(index_x)):
        maxarr.append([index_x[i], index_y[i]])

    print(f'Number of maximum alignments = {len(maxarr)}')
    print('-' * 40)

    for i in range(len(maxarr)):
        max_inx = maxarr[i][0]
        max_iny = maxarr[i][1]
        score = highmax
        firstAlignment = ''
        secondAlignment = ''
        while score > 0:
            if Matrix[max_inx][max_iny] == Matrix[max_inx - 1][max_iny - 1] + match and \
                    S1[max_iny - 1] == S2[max_inx - 1]:
                firstAlignment = S1[max_iny - 1] + firstAlignment
                secondAlignment = S2[max_inx - 1] + secondAlignment
                score = score - match
                max_inx = max_inx - 1
                max_iny = max_iny - 1
            elif Matrix[max_inx][max_iny] == Matrix[max_inx - 1][max_iny - 1] + mismatch:
                firstAlignment = S1[max_iny - 1] + firstAlignment
                secondAlignment = S2[max_inx - 1] + secondAlignment
                score = score - mismatch
                max_inx = max_inx - 1
                max_iny = max_iny - 1
            elif Matrix[max_inx][max_iny] == Matrix[max_inx - 1][max_iny] + gap:
                firstAlignment = '-' + firstAlignment
                secondAlignment = S2[max_inx - i] + secondAlignment
                score = score - gap
                max_inx = max_inx - 1
            elif Matrix[max_inx][max_iny] == Matrix[max_inx][max_iny - 1] + gap:
                firstAlignment = S1[max_iny - 1] + firstAlignment
                secondAlignment = '-' + secondAlignment
                score = score - gap
                max_iny = max_iny - 1

        print(f'Local Alignment')
        print(f"Align1: {firstAlignment}")
        print(f"Align2: {secondAlignment}")
        print('-' * 40)


def globalAlignmentForProtein(X, Y):
    Blosum62 = {"AA": 4, "AC": 0, "AD": -2, "AE": -1, "AF": -2, "AG": 0, "AH": -2, "AI": -1, "AK": -1, "AL": -1,
                "AM": -1,
                "AN": -2, "AP": -1, "AQ": -1, "AR": -1, "AS": 1, "AT": -1, "AV": 0, "AW": -3, "AY": -2, "CA": 0,
                "CC": 9,
                "CD": -3, "CE": -4, "CF": -2, "CG": -3, "CH": -3, "CI": -1, "CK": -3, "CL": -1, "CM": -1, "CN": -3,
                "CP": -3, "CQ": -3, "CR": -3, "CS": -1, "CT": -1, "CV": -1, "CW": -2, "CY": -2,
                "DA": -2, "DC": -3, "DD": 6, "DE": 2, "DF": -3, "DG": -1, "DH": 1, "DI": -3, "DK": -1, "DL": -4,
                "DM": -3,
                "DN": 1, "DP": -1, "DQ": 0, "DR": -2, "DS": 0, "DT": 1, "DV": -3, "DW": -4, "DY": -3,
                "EA": -1, "EC": -4, "ED": 2, "EE": 5, "EF": -3, "EG": -2, "EH": 0, "EI": -3, "EK": 1, "EL": -3,
                "EM": -2,
                "EN": 0, "EP": -1, "EQ": 2, "ER": 0, "ES": 0, "ET": 0, "EV": -2, "EW": -3, "EY": -2,
                "FA": -2, "FC": -2, "FD": -3, "FE": -3, "FF": 6, "FG": -3, "FH": -1, "FI": 0, "FK": -3, "FL": 0,
                "FM": 0,
                "FN": -3, "FP": -4, "FQ": -3, "FR": -3, "FS": -2, "FT": -2, "FV": -1, "FW": 1, "FY": 3,
                "GA": 0, "GC": -3, "GD": -1, "GE": -2, "GF": -3, "GG": 6, "GH": -2, "GI": -4, "GK": -2, "GL": -4,
                "GM": -3,
                "GN": 0, "GP": -2, "GQ": -2, "GR": -2, "GS": 0, "GT": 1, "GV": -3, "GW": -2, "GY": -3,
                "HA": -2, "HC": -3, "HD": -1, "HE": 0, "HF": -1, "HG": -2, "HH": 8, "HI": -3, "HK": -1, "HL": -3,
                "HM": -2,
                "HN": -1, "HP": -2, "HQ": 0, "HR": 0, "HS": -1, "HT": 0, "HV": -3, "HW": -2, "HY": 2,
                "IA": -1, "IC": -1, "ID": -3, "IE": -3, "IF": 0, "IG": -4, "IH": -3, "II": 4, "IK": -3, "IL": 2,
                "IM": 1,
                "IN": -3, "IP": -3, "IQ": -3, "IR": -3, "IS": -2, "IT": -2, "IV": 3, "IW": -3, "IY": -1,
                "KA": -1, "KC": -3, "KD": -1, "KE": 1, "KF": -3, "KG": -2, "KH": -1, "KI": -3, "KK": 5, "KL": -2,
                "KM": -1, "KN": 0, "KP": -1, "KQ": 1, "KR": 2, "KS": 0, "KT": 0, "KV": 0, "KW": -2, "KY": -3,
                "LA": -1, "LC": -1, "LD": -4, "LE": -3, "LF": 0, "LG": -4, "LH": -3, "LI": 2, "LK": -2, "LL": 4,
                "LM": 2,
                "LN": -3, "LP": -3, "LQ": -2, "LR": -2, "LS": -2, "LT": -2, "LV": 1, "LW": -2, "LY": -1,
                "MA": -1, "MC": -1, "MD": -3, "ME": -2, "MF": 0, "MG": -3, "MH": -2, "MI": 1, "MK": -1, "ML": 2,
                "MM": 5,
                "MN": -2, "MP": -2, "MQ": 0, "MR": -1, "MS": -1, "MT": -1, "MV": 1, "MW": -1, "MY": -1,
                "NA": -1, "NC": -3, "ND": 1, "NE": 0, "NF": -3, "NG": -2, "NH": 1, "NI": -3, "NK": 0, "NL": -3,
                "NM": -2,
                "NN": 6, "NP": -1, "NQ": 0, "NR": 0, "NS": 1, "NT": 0, "NV": -3, "NW": -4, "NY": -2,
                "PA": -1, "PC": -3, "PD": -1, "PE": -1, "PF": -4, "PG": -2, "PH": -2, "PI": -3, "PK": -1, "PL": -3,
                "PM": -2,
                "PN": -2, "PP": 7, "PQ": -1, "PR": -2, "PS": -1, "PT": 1, "PV": -2, "PW": -4, "PY": -3,
                "QA": -1, "QC": -3, "QD": 0, "QE": 2, "QF": -3, "QG": -2, "QH": 0, "QI": -3, "QK": 1, "QL": -2, "QM": 0,
                "QN": 0, "QP": -1, "QQ": 5, "QR": 1, "QS": 0, "QT": 0, "QV": -2, "QW": -2, "QY": -1,
                "RA": -1, "RC": -3, "RD": -2, "RE": 0, "RF": -3, "RG": -2, "RH": 0, "RI": -3, "RK": 2, "RL": -2,
                "RM": -1,
                "RN": 0, "RP": -2, "RQ": 1, "RR": 5, "RS": -1, "RT": -1, "RV": -3, "RW": -3, "RY": -2,
                "SA": 1, "SC": -1, "SD": 0, "SE": 0, "SF": -2, "SG": 0, "SH": -1, "SI": -2, "SK": 0, "SL": -2, "SM": -1,
                "SN": 1, "SP": -1, "SQ": 0, "SR": -1, "SS": 4, "ST": 1, "SV": -2, "SW": -3, "SY": -2,
                "TA": -1, "TC": -1, "TD": 1, "TE": 0, "TF": -2, "TG": 1, "TH": 0, "TI": -2, "TK": 0, "TL": -2, "TM": -1,
                "TN": 0, "TP": 1, "TQ": 0, "TR": -1, "TS": 1, "TT": 5, "TV": -2, "TW": -3, "TY": -2,
                "VA": -2, "VC": -1, "VD": -3, "VE": -3, "VF": -1, "VG": 0, "VH": -2, "VI": 1, "VK": -3, "VL": 3,
                "VM": -2,
                "VN": -3, "VP": -2, "VQ": -2, "VR": -3, "VS": -2, "VT": -2, "VV": 4, "VW": -3, "VY": -1,
                "WA": -3, "WC": -2, "WD": -4, "WE": -3, "WF": 1, "WG": -2, "WH": -2, "WI": -3, "WK": -3, "WL": -2,
                "WM": -1,
                "WN": -4, "WP": -4, "WQ": -2, "WR": -3, "WS": -3, "WT": 0, "WV": -3, "WW": 11, "WY": 2,
                "YA": -2, "YC": -2, "YD": -3, "YE": -2, "YF": 3, "YG": -3, "YH": 2, "YI": -1, "YK": -2, "YL": -1,
                "YM": -1,
                "YN": -2, "YP": -3, "YQ": -1, "YR": -2, "YS": -2, "YT": -1, "YV": -1, "YW": 2, "YY": 7}
    m = len(X)
    n = len(Y)
    l = [[0 for X in range(m + 1)] for Y in range(n + 1)]
    for i in range(n + 1):
        l[i][0] = -i
    for j in range(m + 1):
        l[0][j] = -j
    for i in range(1, n + 1, 1):
        for j in range(1, m + 1, 1):
            l[i][j] = max(l[i - 1][j] - 1, l[i][j - 1] - 1, l[i - 1][j - 1] + Blosum62[X[j - 1] + Y[i - 1]])
    for line in l:
        print(line)

    gap1 = ""
    gap2 = ""
    match = ""
    i = n
    j = m
    while i > 0 and j > 0:
        up_left = l[i - 1][j - 1]
        up = l[i - 1][j]

        score = Blosum62[X[j - 1] + Y[i - 1]]
        diagonal = up_left + score
        if l[i][j] == diagonal:
            gap1 += X[j - 1]
            gap2 += Y[i - 1]

            if X[j - 1] == Y[i - 1]:
                match += "|"
            else:
                match += " "
            i = i - 1
            j = j - 1
        elif l[i][j] == up - 1:
            gap1 += "-"
            match += " "
            gap2 += Y[i - 1]
            i = i - 1
        else:
            gap1 += X[j - 1]
            gap2 += "-"
            match += " "
            j = j - 1

    while i > 0:
        gap1 += "-"
        gap2 += Y[i - 1]
        match += " "
        i = i - 1
    while j > 0:
        gap1 += X[j - 1]
        gap2 += "-"
        match += " "
        j = j - 1
    gap1 = gap1[::-1]
    gap2 = gap2[::-1]
    match = match[::-1]

    print(gap1)
    print(match)
    print(gap2)


def localAlignmentForProtein(x, y):
    Blosum62 = {"AA": 4, "AC": 0, "AD": -2, "AE": -1, "AF": -2, "AG": 0, "AH": -2, "AI": -1, "AK": -1, "AL": -1,
                "AM": -1,
                "AN": -2, "AP": -1, "AQ": -1, "AR": -1, "AS": 1, "AT": -1, "AV": 0, "AW": -3, "AY": -2, "CA": 0,
                "CC": 9,
                "CD": -3, "CE": -4, "CF": -2, "CG": -3, "CH": -3, "CI": -1, "CK": -3, "CL": -1, "CM": -1, "CN": -3,
                "CP": -3, "CQ": -3, "CR": -3, "CS": -1, "CT": -1, "CV": -1, "CW": -2, "CY": -2,
                "DA": -2, "DC": -3, "DD": 6, "DE": 2, "DF": -3, "DG": -1, "DH": 1, "DI": -3, "DK": -1, "DL": -4,
                "DM": -3,
                "DN": 1, "DP": -1, "DQ": 0, "DR": -2, "DS": 0, "DT": 1, "DV": -3, "DW": -4, "DY": -3,
                "EA": -1, "EC": -4, "ED": 2, "EE": 5, "EF": -3, "EG": -2, "EH": 0, "EI": -3, "EK": 1, "EL": -3,
                "EM": -2,
                "EN": 0, "EP": -1, "EQ": 2, "ER": 0, "ES": 0, "ET": 0, "EV": -2, "EW": -3, "EY": -2,
                "FA": -2, "FC": -2, "FD": -3, "FE": -3, "FF": 6, "FG": -3, "FH": -1, "FI": 0, "FK": -3, "FL": 0,
                "FM": 0,
                "FN": -3, "FP": -4, "FQ": -3, "FR": -3, "FS": -2, "FT": -2, "FV": -1, "FW": 1, "FY": 3,
                "GA": 0, "GC": -3, "GD": -1, "GE": -2, "GF": -3, "GG": 6, "GH": -2, "GI": -4, "GK": -2, "GL": -4,
                "GM": -3,
                "GN": 0, "GP": -2, "GQ": -2, "GR": -2, "GS": 0, "GT": 1, "GV": -3, "GW": -2, "GY": -3,
                "HA": -2, "HC": -3, "HD": -1, "HE": 0, "HF": -1, "HG": -2, "HH": 8, "HI": -3, "HK": -1, "HL": -3,
                "HM": -2,
                "HN": -1, "HP": -2, "HQ": 0, "HR": 0, "HS": -1, "HT": 0, "HV": -3, "HW": -2, "HY": 2,
                "IA": -1, "IC": -1, "ID": -3, "IE": -3, "IF": 0, "IG": -4, "IH": -3, "II": 4, "IK": -3, "IL": 2,
                "IM": 1,
                "IN": -3, "IP": -3, "IQ": -3, "IR": -3, "IS": -2, "IT": -2, "IV": 3, "IW": -3, "IY": -1,
                "KA": -1, "KC": -3, "KD": -1, "KE": 1, "KF": -3, "KG": -2, "KH": -1, "KI": -3, "KK": 5, "KL": -2,
                "KM": -1, "KN": 0, "KP": -1, "KQ": 1, "KR": 2, "KS": 0, "KT": 0, "KV": 0, "KW": -2, "KY": -3,
                "LA": -1, "LC": -1, "LD": -4, "LE": -3, "LF": 0, "LG": -4, "LH": -3, "LI": 2, "LK": -2, "LL": 4,
                "LM": 2,
                "LN": -3, "LP": -3, "LQ": -2, "LR": -2, "LS": -2, "LT": -2, "LV": 1, "LW": -2, "LY": -1,
                "MA": -1, "MC": -1, "MD": -3, "ME": -2, "MF": 0, "MG": -3, "MH": -2, "MI": 1, "MK": -1, "ML": 2,
                "MM": 5,
                "MN": -2, "MP": -2, "MQ": 0, "MR": -1, "MS": -1, "MT": -1, "MV": 1, "MW": -1, "MY": -1,
                "NA": -1, "NC": -3, "ND": 1, "NE": 0, "NF": -3, "NG": -2, "NH": 1, "NI": -3, "NK": 0, "NL": -3,
                "NM": -2,
                "NN": 6, "NP": -1, "NQ": 0, "NR": 0, "NS": 1, "NT": 0, "NV": -3, "NW": -4, "NY": -2,
                "PA": -1, "PC": -3, "PD": -1, "PE": -1, "PF": -4, "PG": -2, "PH": -2, "PI": -3, "PK": -1, "PL": -3,
                "PM": -2,
                "PN": -2, "PP": 7, "PQ": -1, "PR": -2, "PS": -1, "PT": 1, "PV": -2, "PW": -4, "PY": -3,
                "QA": -1, "QC": -3, "QD": 0, "QE": 2, "QF": -3, "QG": -2, "QH": 0, "QI": -3, "QK": 1, "QL": -2, "QM": 0,
                "QN": 0, "QP": -1, "QQ": 5, "QR": 1, "QS": 0, "QT": 0, "QV": -2, "QW": -2, "QY": -1,
                "RA": -1, "RC": -3, "RD": -2, "RE": 0, "RF": -3, "RG": -2, "RH": 0, "RI": -3, "RK": 2, "RL": -2,
                "RM": -1,
                "RN": 0, "RP": -2, "RQ": 1, "RR": 5, "RS": -1, "RT": -1, "RV": -3, "RW": -3, "RY": -2,
                "SA": 1, "SC": -1, "SD": 0, "SE": 0, "SF": -2, "SG": 0, "SH": -1, "SI": -2, "SK": 0, "SL": -2, "SM": -1,
                "SN": 1, "SP": -1, "SQ": 0, "SR": -1, "SS": 4, "ST": 1, "SV": -2, "SW": -3, "SY": -2,
                "TA": -1, "TC": -1, "TD": 1, "TE": 0, "TF": -2, "TG": 1, "TH": 0, "TI": -2, "TK": 0, "TL": -2, "TM": -1,
                "TN": 0, "TP": 1, "TQ": 0, "TR": -1, "TS": 1, "TT": 5, "TV": -2, "TW": -3, "TY": -2,
                "VA": -2, "VC": -1, "VD": -3, "VE": -3, "VF": -1, "VG": 0, "VH": -2, "VI": 1, "VK": -3, "VL": 3,
                "VM": -2,
                "VN": -3, "VP": -2, "VQ": -2, "VR": -3, "VS": -2, "VT": -2, "VV": 4, "VW": -3, "VY": -1,
                "WA": -3, "WC": -2, "WD": -4, "WE": -3, "WF": 1, "WG": -2, "WH": -2, "WI": -3, "WK": -3, "WL": -2,
                "WM": -1,
                "WN": -4, "WP": -4, "WQ": -2, "WR": -3, "WS": -3, "WT": 0, "WV": -3, "WW": 11, "WY": 2,
                "YA": -2, "YC": -2, "YD": -3, "YE": -2, "YF": 3, "YG": -3, "YH": 2, "YI": -1, "YK": -2, "YL": -1,
                "YM": -1,
                "YN": -2, "YP": -3, "YQ": -1, "YR": -2, "YS": -2, "YT": -1, "YV": -1, "YW": 2, "YY": 7}
    m = len(x)
    n = len(y)
    L = [[0 for x in range(m + 1)] for y in range(n + 1)]
    m1 = 0
    m2 = 0
    maxNumber = -1
    for i in range(n + 1):
        L[i][0] = 0
    for j in range(m + 1):
        L[0][j] = 0
    for i in range(1, n + 1, 1):
        for j in range(1, m + 1, 1):

            L[i][j] = max(L[i - 1][j] - 1, 0, L[i][j - 1] - 1, L[i - 1][j - 1] + Blosum62[x[j - 1] + y[i - 1]])
            if (L[i][j] > maxNumber):
                maxNumber = L[i][j]
                m1 = i
                m2 = j
    for line in L:
        print(line)
    # print(m1, m2)
    # print(L[m1][m2])

    String1 = ""
    String2 = ""
    Match = ""
    i = m1
    j = m2
    while i > 0 and j > 0:
        up_left = L[i - 1][j - 1]
        up = L[i - 1][j]
        score = Blosum62[x[j - 1] + y[i - 1]]
        diagonal = up_left + score
        if L[i][j] == diagonal:
            String1 += x[j - 1]
            String2 += y[i - 1]
            if x[j - 1] == y[i - 1]:
                Match += "|"
            else:
                Match += " "
            i = i - 1
            j = j - 1
        elif L[i][j] == up - 1:
            String1 += "-"
            Match += " "
            String2 += y[i - 1]

            i = i - 1
        else:
            String1 += x[j - 1]
            String2 += "-"
            Match += " "
            j = j - 1
    while (i > 0):
        String1 += "-"
        String2 += y[i - 1]
        Match += " "
        i = i - 1
    while (j > 0):
        String1 += x[j - 1]
        String2 += "-"
        Match += " "
        j = j - 1

    String1 = String1[::-1]
    String2 = String2[::-1]
    Match = Match[::-1]

    print(String1)
    print(Match)
    print(String2)


if __name__ == '__main__':
    print("Choose Type Alignment For:")
    print("1- Global Alignment For DNA")
    print("2- Local Alignment For DNA")
    print("3- Global Alignment For Protein")
    print("4- Local Alignment For Protein")
    print('*' * 40)
    choose = int(input("Enter Your Choice: "))
    if choose == 1:
        globalAlignmentForDNA()
    elif choose == 2:
        localAlignmentForDNA()
    elif choose == 3:
        Seq1 = input("Enter Sequence 1: ").upper()
        Seg2 = input("Enter Sequence 2: ").upper()
        globalAlignmentForProtein(Seq1, Seg2)
    else:
        Seq1 = input("Enter Sequence 1: ").upper()
        Seg2 = input("Enter Sequence 2: ").upper()
        localAlignmentForProtein(Seq1, Seg2)
