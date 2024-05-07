def nucl_complement(sequence):
    trantab = str.maketrans("ACGTacgt", "TGCAtgca")
    transtr = sequence.translate(trantab)
    return transtr[::-1]


def get_coding_seq(dna_seq):
    start_codons = ["ATG", "TTG", "GTG"]
    for i in range(len(dna_seq) - 2):
        if dna_seq[i : i + 3] in start_codons:
            coding_seq = dna_seq[i:]
            if len(coding_seq) % 3 == 0:
                return coding_seq
            else:
                return coding_seq[: len(coding_seq) // 3 * 3]
