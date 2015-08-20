from converter_fasta import convert_fasta_to_string
from genome_classes import Genome, Domain, Target


def init_genome_camv(dt=0.001, complex_concentration=135000000000):
    # initialize camv genome
    sequence = convert_fasta_to_string("genome_camv.fasta")
    genome_camv = Genome(sequence)

    # initialize camv genome domains
    #ncr_left_of_35S = ...
    #ncr_right_of_35S = ...
    #promoter_35S = ...
    #promoter_19S = ...
    #gene_P1 = ...
    #gene_P2 = ...
    #gene_P3 = ...
    #gene_P4 = ...
    #gene_P5 = ...
    gene_P6 = Domain("P6", 500, 600, "orf", genome_camv)

    # initialize camv genome targets
    target_P6_1 = Target("P6_1", "ggagaaagaaaagatatttaaaa", "ggagaaagaaaagatatttaaaa", 521, complex_concentration, 1, dt, gene_P6)

    """ isaac's testing block
    print P6_1.grna
    print P6.target_location("P6_1")
    print P6_1.cut_probability

    P6_1.cut()
    P6_1.repair(dt)

    print P6_1.sequence
    print P6.target_location("P6_1")
    print P6_1.cut_probability
    print P6_1.shift
    """

    return genome_camv
