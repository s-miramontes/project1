# DNA -> RNA Transcription


def transcribe(seq: str) -> str:
    """
    transcribes DNA to RNA by replacing
    all `T` to `U`
    """

    # rna dict
    rna_dict = {'A': 'U', 'T': 'A',
                'C': 'G', 'G': 'C'}

    # seq to list
    seq_list = list(seq)

    # for each el in list get key
    for i, s in enumerate(seq_list):
        seq_list[i] = rna_dict[s]

    # list to string
    rna_seq = ''.join(seq_list)
    
    return rna_seq


def reverse_transcribe(seq: str) -> str:
    """
    transcribes DNA to RNA by replacing
    all `T` to `U` then reverses the sequence
    """
    # Change T to U first
    transc = transcribe(seq)

    # return the reverse order
    return transc[::-1]


