import parasail as ps

MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGT", "TGCA")


def reverse_and_complement(query: str) -> str:
    """Makes a complement and reverses the given query sequence.
    :param query: A query DNA sequence.
    """
    return query.translate(COMP)[::-1]


# TODO: In essence this has turned into a largest common substring problem instead of an alignment problem. What if I allow at most math.ceil(len(sequence)/10) gaps?
def align(query: str, seq: str, gap_open: int = 1) -> dict[str, int | str]:
    """Performs Smith-Waterman on the given query and sequence.
    :param query: DNA sequence 1.
    :param seq: DNA sequence 2.
    :param gap_open: penalty for opening and extending a gap.
    """
    result = ps.sw_trace(query, seq, gap_open, gap_open, MATRIX)
    aligned_primer = result.traceback.ref
    aligned_primer = aligned_primer.replace("-", "N")
    return {"alignment score": result.score, "aligned portion": aligned_primer}
