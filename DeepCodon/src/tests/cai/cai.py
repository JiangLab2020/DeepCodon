def CAIC(sequence, condonAA, RSCUNN, RSCUNN1):
    CAI0 = 1
    n = 0
    for i in range(0, len(sequence), 3):
        if i + 3 < len(sequence):
            sequencesTTT = sequence[i : i + 3].upper()
            CAI0 = CAI0 * (RSCUNN1[sequencesTTT] / RSCUNN[condonAA[sequencesTTT]])
            n += 1
    CAIN = CAI0 ** (1 / n)
    return CAIN


def processing(infile, inseq):
    with open(infile) as f:
        records = f.read()
    records = records.split()
    encodings = []
    encodings0 = []
    encodings1 = []
    encodings2 = []
    for i in range(0, 320, 5):
        encodings.append(records[i])
        encodings.append(records[i + 1])
        encodings.append(records[i + 2])
        encodings0.append(records[i])
        encodings1.append(records[i + 1])
        encodings2.append(int(float(records[i + 2])))
    MM = [
        "Gly",
        "Glu",
        "Asp",
        "Val",
        "Ala",
        "Arg",
        "Ser",
        "Lys",
        "Asn",
        "Met",
        "Ile",
        "Thr",
        "Trp",
        "End",
        "Cys",
        "Tyr",
        "Leu",
        "Phe",
        "Gln",
        "His",
        "Pro",
    ]
    values1 = [0] * 21
    values2 = ["QQQ"] * 21
    values3 = [0] * 64

    RSCUFM = dict(zip(MM, values1))

    RSCUFZ = dict(zip(MM, values1))

    RSCUFZ1 = dict(zip(encodings1, encodings2))

    RSCUFZA = dict(zip(MM, values2))

    RSCUNN = dict(zip(MM, values1))

    RSCUNN1 = dict(zip(encodings1, values3))

    RSCUFZ_LIST = []
    for i in range(42):
        RSCUFZ_LIST.append([])
    for i in range(0, 192):
        k = 0
        for j in MM:
            if encodings[i] == j:
                str1 = str(encodings[i + 1])
                m = int(float(encodings[i + 2]))
                RSCUFM[j] = RSCUFM[j] + m
                RSCUFZ_LIST[k].append(str1)
                RSCUFZ_LIST[k + 1].append(m)
            k += 2

    mmm = [4, 2, 2, 4, 4, 6, 6, 2, 2, 1, 3, 4, 1, 3, 2, 2, 6, 2, 2, 2, 4]

    j = 0
    k = 0
    for i in MM:
        n = RSCUFZ_LIST[j + 1].index(max(RSCUFZ_LIST[j + 1]))
        RSCUFZA[i] = str(RSCUFZ_LIST[j][n])
        RSCUFZ[i] = int(RSCUFZ_LIST[j + 1][n]) * mmm[k]
        j = j + 2
        k = k + 1
    for i in MM:
        RSCUNN[i] = int(RSCUFZ[i]) / int(RSCUFM[i])

    AA = "GATC"
    triN = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]

    mmm1 = [
        0,
        0,
        0,
        0,
        1,
        1,
        2,
        2,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        5,
        5,
        6,
        6,
        7,
        7,
        8,
        8,
        9,
        10,
        10,
        10,
        11,
        11,
        11,
        11,
        12,
        13,
        14,
        14,
        13,
        13,
        15,
        15,
        16,
        16,
        17,
        17,
        6,
        6,
        6,
        6,
        5,
        5,
        5,
        5,
        18,
        18,
        19,
        19,
        16,
        16,
        16,
        16,
        20,
        20,
        20,
        20,
    ]

    mmm2 = [
        4,
        4,
        4,
        4,
        2,
        2,
        2,
        2,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        6,
        6,
        6,
        6,
        2,
        2,
        2,
        2,
        1,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        1,
        3,
        2,
        2,
        3,
        3,
        2,
        2,
        6,
        6,
        2,
        2,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        6,
        2,
        2,
        2,
        2,
        6,
        6,
        6,
        6,
        4,
        4,
        4,
        4,
    ]

    for i in range(64):
        RSCUNN1[triN[i]] = (RSCUFZ1[triN[i]] / RSCUFM[MM[mmm1[i]]]) * int(mmm2[i])

    condonAA = dict(zip(encodings1, encodings0))

    caiout = CAIC(inseq, condonAA, RSCUNN, RSCUNN1)
    return caiout
