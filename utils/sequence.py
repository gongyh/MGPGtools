def nucl_complement(sequence):
    trantab = str.maketrans("ACGTacgt", "TGCAtgca")
    transtr = sequence.translate(trantab)
    return transtr[::-1]


def get_protein_coding_regions(dna_seq):
    start_codons = ["ATG", "TTG", "TGT"]
    stop_codons = ["TAA", "TAG", "TGA"]
    protein_coding_seq = ""
    started = False
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i : i + 3]
        if codon in start_codons:
            started = True
        if started:
            protein_coding_seq += codon
        if codon in stop_codons and started:
            break
    if len(protein_coding_seq) % 3 == 0:
        return protein_coding_seq
    else:
        return protein_coding_seq[: len(protein_coding_seq) // 3 * 3]
