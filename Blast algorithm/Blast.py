import os
#Ahmed Tarek Elkomy
#Mostafa mohamed mahmoud
#Marwan Akram Said
#Moahmed Afifi Abd el hafez
#MMohamed Ahmed Moahmed Moastafa
# Step 1
user_query = input('Enter a query:- ').upper()
sep = [user_query[i:i + 2] for i in range(0, len(user_query), 2)]
path = os.path.dirname(os.path.abspath(__file__))
File = open(path + "/Database.txt", "r")
Database = [line.rstrip() for line in File.readlines()]
File.close()
Blosum62 = {"AA": 4, "AC": 0, "AD": -2, "AE": -1, "AF": -2, "AG": 0, "AH": -2, "AI": -1, "AK": -1, "AL": -1, "AM": -1,
            "AN": -2, "AP": -1, "AQ": -1, "AR": -1, "AS": 1, "AT": -1, "AV": 0, "AW": -3, "AY": -2, "CA": 0, "CC": 9,
            "CD": -3, "CE": -4, "CF": -2, "CG": -3, "CH": -3, "CI": -1, "CK": -3, "CL": -1, "CM": -1, "CN": -3,
            "CP": -3, "CQ": -3, "CR": -3, "CS": -1, "CT": -1, "CV": -1, "CW": -2, "CY": -2,
            "DA": -2, "DC": -3, "DD": 6, "DE": 2, "DF": -3, "DG": -1, "DH": 1, "DI": -3, "DK": -1, "DL": -4, "DM": -3,
            "DN": 1, "DP": -1, "DQ": 0, "DR": -2, "DS": 0, "DT": 1, "DV": -3, "DW": -4, "DY": -3,
            "EA": -1, "EC": -4, "ED": 2, "EE": 5, "EF": -3, "EG": -2, "EH": 0, "EI": -3, "EK": 1, "EL": -3, "EM": -2,
            "EN": 0, "EP": -1, "EQ": 2, "ER": 0, "ES": 0, "ET": -1, "EV": -2, "EW": -3, "EY": -2,
            "FA": -2, "FC": -2, "FD": -3, "FE": -3, "FF": 6, "FG": -3, "FH": -1, "FI": 0, "FK": -3, "FL": 0, "FM": 0,
            "FN": -3, "FP": -4, "FQ": -3, "FR": -3, "FS": -2, "FT": -2, "FV": -1, "FW": 1, "FY": 3,
            "GA": 0, "GC": -3, "GD": -1, "GE": -2, "GF": -3, "GG": 6, "GH": -2, "GI": -4, "GK": -2, "GL": -4,
            "GM": -3,
            "GN": 0, "GP": -2, "GQ": -2, "GR": -2, "GS": 0, "GT": 1, "GV": -3, "GW": -2, "GY": -3,
            "HA": -2, "HC": -3, "HD": -1, "HE": 0, "HF": -1, "HG": -2, "HH": 8, "HI": -3, "HK": -1, "HL": -3,
            "HM": -2,
            "HN": -1, "HP": -2, "HQ": 0, "HR": 0, "HS": -1, "HT": 0, "HV": -3, "HW": -2, "HY": 2,
            "IA": -1, "IC": -1, "ID": -3, "IE": -3, "IF": 0, "IG": -4, "IH": -3, "II": 4, "IK": -3, "IL": 2, "IM": 1,
            "IN": -3, "IP": -3, "IQ": -3, "IR": -3, "IS": -2, "IT": -2, "IV": 3, "IW": -3, "IY": -1,
            "KA": -1, "KC": -3, "KD": -1, "KE": 1, "KF": -3, "KG": -2, "KH": -1, "KI": -3, "KK": 5, "KL": -2,
            "KM": -1, "KN": 0, "KP": -1, "KQ": 1, "KR": 2, "KS": 0, "KT": 0, "KV": 0, "KW": -2, "KY": -3,
            "LA": -1, "LC": -1, "LD": -4, "LE": -3, "LF": 0, "LG": -4, "LH": -3, "LI": 2, "LK": -2, "LL": 4, "LM": 2,
            "LN": -3, "LP": -3, "LQ": -2, "LR": -2, "LS": -2, "LT": -2, "LV": 1, "LW": -2, "LY": -1,
            "MA": -1, "MC": -1, "MD": -3, "ME": -2, "MF": 0, "MG": -3, "MH": -2, "MI": 1, "MK": -1, "ML": 2, "MM": 5,
            "MN": -2, "MP": -2, "MQ": 0, "MR": -1, "MS": -1, "MT": -1, "MV": 1, "MW": -1, "MY": -1,
            "NA": -1, "NC": -3, "ND": 1, "NE": 0, "NF": -3, "NG": -2, "NH": 1, "NI": -3, "NK": 0, "NL": -3, "NM": -2,
            "NN": 6, "NP": -1, "NQ": 0, "NR": 0, "NS": 1, "NT": 0, "NV": -3, "NW": -4, "NY": -2,
            "PA": -1, "PC": -3, "PD": -1, "PE": -1, "PF": -4, "PG": -2, "PH": -2, "PI": -3, "PK": -1, "PL": -3,
            "PM": -2,
            "PN": -2, "PP": 7, "PQ": -1, "PR": -2, "PS": -1, "PT": 1, "PV": -2, "PW": -4, "PY": -3,
            "QA": -1, "QC": -3, "QD": 0, "QE": 2, "QF": -3, "QG": -2, "QH": 0, "QI": -3, "QK": 1, "QL": -2, "QM": 0,
            "QN": 0, "QP": -1, "QQ": 5, "QR": 1, "QS": 0, "QT": 0, "QV": -2, "QW": -2, "QY": -1,
            "RA": -1, "RC": -3, "RD": -2, "RE": 0, "RF": -3, "RG": -2, "RH": 0, "RI": -3, "RK": 2, "RL": -2, "RM": -1,
            "RN": 0, "RP": -2, "RQ": 1, "RR": 5, "RS": -1, "RT": -1, "RV": -3, "RW": -3, "RY": -2,
            "SA": 1, "SC": -1, "SD": 0, "SE": 0, "SF": -2, "SG": 0, "SH": -1, "SI": -2, "SK": 0, "SL": -2, "SM": -1,
            "SN": 1, "SP": -1, "SQ": 0, "SR": -1, "SS": 4, "ST": 1, "SV": -2, "SW": -3, "SY": -2,
            "TA": -1, "TC": -1, "TD": 1, "TE": -1, "TF": -2, "TG": 1, "TH": 0, "TI": -2, "TK": 0, "TL": -2, "TM": -1,
            "TN": 0, "TP": 1, "TQ": 0, "TR": -1, "TS": 1, "TT": 5, "TV": -2, "TW": -3, "TY": -2,
            "VA": -2, "VC": -1, "VD": -3, "VE": -3, "VF": -1, "VG": 0, "VH": -2, "VI": 1, "VK": -3, "VL": 3, "VM": -2,
            "VN": -3, "VP": -2, "VQ": -2, "VR": -3, "VS": -2, "VT": -2, "VV": 4, "VW": -3, "VY": -1,
            "WA": -3, "WC": -2, "WD": -4, "WE": -3, "WF": 1, "WG": -2, "WH": -2, "WI": -3, "WK": -3, "WL": -2,
            "WM": -1,
            "WN": -4, "WP": -4, "WQ": -2, "WR": -3, "WS": -3, "WT": 0, "WV": -3, "WW": 11, "WY": 2,
            "YA": -2, "YC": -2, "YD": -3, "YE": -2, "YF": 3, "YG": -3, "YH": 2, "YI": -1, "YK": -2, "YL": -1,
            "YM": -1,
            "YN": -2, "YP": -3, "YQ": -1, "YR": -2, "YS": -2, "YT": -1, "YV": -1, "YW": 2, "YY": 7
    , "AX": -1, "RX": -1, "NX": -1, "DX": -1, "CX": -1, "QX": -1, "EX": -1, "GX": -1, "HX": -1, "IX": -1, "LX": -1,
            "KX": -1, "MX": -1, "FX": -1, "PX": -1, "SX": -1, "TX": -1, "WX": -1, "YX": -1, "VX": -1, "XX": -1, }
# STEP 1 #remove low comlexity
for i in range(0, len(sep) - 1):
    if sep[i] == sep[i + 1]:
        user_query = user_query.replace(sep[i], 'XX')
# STEP 2
wordL = int(input('Enter word lenght :-'))
num = len(user_query) - wordL
wordnum = num + 1
k = []
for i in range(0, wordnum):
    k.append(user_query[i:i + wordL])
# PQGEFG
# STEP 3,4
wthreshold = int(input("enter word threshold :-"))
seeds = dict()
s = "ARNDCQEGHILKMFPSTWYV"

for i in range(0, wordnum):  # 3adad el words
    for j in range(0, wordL):  # lenght of word
        for l in range(0, 20):  # change to other aa
            compare = k[i]
            compare = compare.replace(compare[j], s[l])
            score = 0
            for m in range(0, wordL):
                y = compare[m] + k[i][m]
                score = score + Blosum62[y]
            if score > wthreshold:
                seeds[compare] = k[i]

# PQGEFG\
HSPW = dict()
HSP = int(input("HSP threshhold :- "))
X = 3
#Step 5,6,7
for i in Database:
    for m in range(0, len(i) - wordL):
        makeExtend = i[m:m + wordL]
        for item in seeds:
            for h in range(0, len(user_query) - wordL):
                extend = user_query[h:h + wordL]
                if makeExtend == item and seeds[item] == extend: #Hit in Database string
                    score = 0
                    for n in range(0, wordL):
                        Extend = makeExtend[n] + extend[n]
                        score = score + Blosum62[Extend]
                        max = score
                    if h == 0 or m == 0: #extend Right only
                        Flagout = False
                        FlagX = False
                        DI = m + wordL - 1
                        QI = h + wordL - 1
                        Oldscore = score
                        while True:
                            DI += 1
                            QI += 1
                            Oldscore = score
                            if QI >= len(user_query) or DI >= len(i):
                                Flagout = True
                                break
                            else:
                                ahmed3 = user_query[QI] + i[DI]
                                score += Blosum62[ahmed3]
                            if score > max:
                                max = score
                            else:
                                HSPP = max - score
                                if HSPP > X:
                                    FlagX = True
                                    break
                        if Oldscore > HSP:
                            if Flagout:
                                HSPW[i[m:DI]] = score
                            elif FlagX:
                                HSPW[i[m:DI]] = Oldscore
                    elif h == len(user_query) or m == len(i):#extends right only
                        DI = m
                        QI = h
                        Flagout = False
                        FlagX = False
                        Oldscore = score
                        while True:
                            DI -= 1
                            QI -= 1
                            Oldscore = score
                            if DI < 0 or QI < 0:
                                Flagout = True
                                break
                            else:
                                query1 = user_query[QI] + i[DI]
                                score += Blosum62[query1]
                            if score > max:
                                max = score
                            else:
                                HSPP = max - score
                                if HSPP > X:
                                    FlagX = True
                        if score > HSP:
                            if Flagout:
                                HSPW[i[DI:m + wordL]] = score
                            elif FlagX:
                                HSPW[i[DI:m + wordL]] = Oldscore
                    else:#extends on both dirction
                        Flagout = False
                        FlagX = False
                        rightD = m + wordL - 1
                        rightQ = h + wordL - 1
                        leftD = m
                        leftQ = h
                        while True:
                            rightD += 1
                            rightQ += 1
                            Oldscore = score
                            if rightD >= len(i) or rightQ >= len(user_query):
                                Flagout = True
                                break
                            else:
                                Q1 = user_query[rightQ] + i[rightD]
                                score += Blosum62[Q1]
                            leftD -= 1
                            leftQ -= 1
                            if leftD < 0 or leftQ < 0:
                                Flagout = True
                                break
                            else:
                                Q2 = user_query[leftQ] + i[leftD]
                                score += Blosum62[Q2]
                            if max < score:
                                max = score
                            else:
                                HSPP = max - score
                                if HSPP > X:
                                    FlagX = True
                                    break
                        if Oldscore > HSP:
                            if Flagout:
                                HSPW[i[leftD:rightD]] = score
                            elif FlagX:
                                HSPW[i[leftD + 1:rightD]] = Oldscore

print(HSPW)
file = open(path + "/Database.txt", "r")
flag = 0
index = 0
Keys = list()
for key in HSPW.keys():
    Keys.append(key)
print()
correctWord = None
for line in file:
    if len(Keys) == 1:
        if Keys[0] in line:
            flag = 1
            correctWord = Keys[0]
        else:
            flag = 0
    else:
        for word in range(len(Keys) - 1):
            if Keys[word] in line and Keys[word + 1] in line:
                flag = 1
                correctWord = [Keys[word], Keys[word + 1]]
            elif Keys[word] in line:
                flag = 1
                correctWord = Keys[word]
            elif Keys[word] not in line and Keys[word + 1] in line:
                flag = 1
                correctWord = Keys[word + 1]
    if flag == 0:
        if correctWord is None:
            print()
    else:
        print('String', correctWord, 'Found In Line In Database Id = ', index)
    index = index + 1

file.close()
