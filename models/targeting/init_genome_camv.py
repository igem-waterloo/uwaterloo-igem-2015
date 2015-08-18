from converter_fasta import convert_fasta_to_string
from genome_classes import Genome, Domain, Target


def init_genome_camv():
    # initialize camv genome, domains, targets
    dt = 0.001
    complex_concentration = 135000000000

    sequence = convert_fasta_to_string("genome_camv.fasta")
    genome_camv = Genome(sequence)
    # These next two aren't at all accurate yet
    P6 = Domain("P6", 500, 600, "promoter", genome_camv)
    P6_1 = Target("P6_1", "ggagaaagaaaagatatttaaaa", "ggagaaagaaaagatatttaaaa", 521, complex_concentration, 1, dt, P6)

    print P6_1.grna
    print P6.target_location("P6_1")
    print P6_1.cut_probability

    P6_1.cut()
    P6_1.repair(dt)

    print P6_1.sequence
    print P6.target_location("P6_1")
    print P6_1.cut_probability
    print P6_1.shift
    return genome_camv
