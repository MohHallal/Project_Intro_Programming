# DNA Translator

## Description

The DNA Translator script is a tool that allows you to convert DNA sequences into RNA sequences and/or amino acid sequences. This translator is coded in Python and provides the following features:

- Conversion between different levels: DNA to RNA, RNA to protein, DNA to protein, based on user choice.
- User can input one or multiple sequences for transcription/translation or provide a FASTA file containing multiple sequences (the FASTA file is validated for its extension, line length, and each base).
- User can also provide a FASTA file with a long nucleic sequence (e.g., a chromosome) accompanied by a GTF/GFF file containing the positions of the genes of interest to be transcribed and/or translated.
- If the expected sequence format (e.g., DNA or RNA) does not match the provided sequence (e.g., RNA), an error will be reported.
- Translation stops when a stop codon is encountered.
- User can provide an alternative codon usage table instead of the standard one (e.g., using Arginine for the codon AUG instead of Methionine).
