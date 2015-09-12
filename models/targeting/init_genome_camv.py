from converter_fasta import convert_fasta_to_string
from genome_classes import Genome, Domain, Target


def init_genome_camv(pseudo_targets):
    """ Initializes the CaMV genome object for simulation
    Args:
        pseudo_targets: list of format [t1, t2, ..., tn]
                        where each is a dict with = {'label': ,
                                                     'target_sequence': ,
                                                     'target_start_idx': ,
                                                     'complex_concentration': ,
                                                     'sense': ,
                                                     'direction': ,
                                                     'domain_label': }
    Returns:
        complete Genome instance corresponding to CaMV with targets
    """
    # initialize camv genome
    sequence = convert_fasta_to_string("genome_camv.fasta")
    genome_camv = Genome(sequence)

    # initialize camv genome 
    # gene coordinates from http://www.ncbi.nlm.nih.gov/nuccore/9626938
    #ncr_left_of_35S = ...
    #genome_camv.add_domain(ncr_left_of_35S)
    #ncr_right_of_35S = ...
    #genome_camv.add_domain(ncr_right_of_35S)
    promoter_35S = Domain("promoter_35S", 7092, 7435, "promoter", genome_camv)
    genome_camv.add_domain(promoter_35S)
    promoter_19S = Domain("promoter_19S", 5380, 5773, "promoter", genome_camv)
    genome_camv.add_domain(promoter_19S)
    gene_P1 = Domain("gene_P1", 364, 1347, "orf", genome_camv)
    genome_camv.add_domain(gene_P1)
    gene_P2 = Domain("gene_P2", 1349, 1828, "orf", genome_camv)
    genome_camv.add_domain(gene_P2)
    gene_P3 = Domain("gene_P3", 1830, 2219, "orf", genome_camv)
    genome_camv.add_domain(gene_P3)
    gene_P4 = Domain("gene_P4", 2201, 3670, "orf", genome_camv)
    genome_camv.add_domain(gene_P4)
    gene_P5 = Domain("gene_P5", 3633, 5672, "orf", genome_camv)
    genome_camv.add_domain(gene_P5)
    gene_P6 = Domain("gene_P6", 5776, 7338, "orf", genome_camv)
    genome_camv.add_domain(gene_P6)

    # initialize camv genome targets
    domain_labels = genome_camv.domains.keys()
    for pseudo_target in pseudo_targets:
        assert pseudo_target.domain_label in domain_labels
        for domain_label in domain_labels:
            if pseudo_target.domain_label == domain_label:
                domain = genome_camv.domains[domain_label]
                target = Target(pseudo_target['label'],
                                pseudo_target['target_sequence'],
                                pseudo_target['target_sequence'],
                                pseudo_target['target_start_idx'],
                                pseudo_target['complex_concentration'],
                                pseudo_target['sense'],
                                pseudo_target['direction'],
                                domain)
                domain.add_target(target)

    """ isaac's testing block
    target_P6_1 = Target("P6_1", "ggagaaagaaaagatatttaaaa", "ggagaaagaaaagatatttaaaa", 521, complex_concentration, 1, 1, gene_P6)
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


def init_targets_single_P6(complex_concentration):
    return [{'label': "target_P6_1",
             'target_sequence': "ggagaaagaaaagatatttaaaa",
             'target_start_idx': 521,
             'complex_concentration': complex_concentration,
             'sense': 1,
             'direction': 1,
             'domain_label': "gene_P6"}]
