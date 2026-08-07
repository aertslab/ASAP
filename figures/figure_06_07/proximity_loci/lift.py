from pyliftover import LiftOver
import pysam


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def liftover_variant(chr:str, bp: int, ref:str, alt:str, fasta, lo):

    REF = ref.upper()
    ALT = alt.upper()

    result = lo.convert_coordinate(chr, bp - 1)

    if not result:
        print(chr, bp, 'failed')
        global fail
        fail = fail + 1
        return 000, REF, ALT

    new_chrom, new_pos0, strand, _ = result[0]
    new_pos = new_pos0 + 1  # back to 1-based

    # If the liftover flips strand, REF/ALT need to be reverse-complemented
    # before comparing against the target's forward-strand sequence.
    exp_ref, exp_alt = (REF, ALT) if strand == "+" else (revcomp(REF), revcomp(ALT))

    fetched_ref = fasta.fetch(new_chrom, new_pos0, new_pos0 + len(exp_ref)).upper()
    
    if fetched_ref == exp_ref:
        return new_pos, exp_ref, exp_alt
    else:
        return new_pos, exp_alt, exp_ref  # swap so new_ref matches the target
