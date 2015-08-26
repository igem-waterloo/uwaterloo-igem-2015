from probabilistic import prob_cut, nt_rand, indel


# Script to store principle simulation classes and their interactions
# Gene/genome de-activation defined by:
#   - frameshift mutation in an ORF
#   - deletion of a promoter
#   - deletion of a significant portion of the gene

# TODO
# - function for cut probability: foo_cut_prob(self.grna, self.target, dt, complex_concentration)
# - function for cut location: foo_cut_posn()
# - function for indel type (L/R del, insert) foo_indel(), foo_nt_rand(insert) etc
# - if repair is big deletion, what do we fill in at the end of the target? cant be random, must come from sequence data
# - code for deletions
# - see specific TODOs / 'to be implemented' throughout the code
# - properly update target.current_start and domain genome repair method


class Target(object):

    def __init__(self, label, grna, sequence, start, complex_concentration, sense, direction, domain):
        self.label = label  # string
        self.grna = grna  # string
        self.sequence = sequence  # string, exclude PAM, should be ~ 20 chars
        self.original_start = start  # int, location of first target nucleotide adjacent to PAM, shouldn't change
        self.current_start = start  # int, location of first target nucleotide adjacent to PAM, changes with indels
        self.total_cuts = 0  # int, total time this target has been cut
        self.cut_position = None  # absolute genome location of cut
        self.repaired = True  # defined by open/closed
        self.targetable = True  # defined by targetable or not (PAM broken or indel size > 5)
        self.complex_concentration = complex_concentration  # conc of gRNA-cas9 complex inside nucleus
        self.shift = 0  # defined by sum of net indel sizes, used to compute frameshift if orf region
        assert sense in [1, -1]
        assert direction in [1, -1]
        self.sense = sense  # 1 or -1, referring to top (explicit) or bottom (implicit) dna strand
        self.direction = direction  # 1 or -1, 1 means ggn occurs before target (numerically)
        self.domain = domain
        domain.add_target(self)
        self.cut_probability = None

    def is_repaired(self):
        return self.repaired

    def is_targetable(self):
        return self.targetable

    def get_shift(self):
        return self.shift

    def compute_and_assign_cut_probability(self, dt):
        self.cut_probability = prob_cut(self.grna, self.sequence, self.complex_concentration, dt)
        return self.cut_probability

    def cut(self):
        self.total_cuts += 1
        self.cut_position = self.current_start + self.direction * 3  # posn of the nt right of the cut, usually 3-4nt from pam: foo_cut_posn()
        self.repaired = False

    def repair(self, dt):
        # call Genome.target_repair through Domain
        net_indel_size, self.sequence = self.domain.genome_repair(self.label, self.cut_position)

        # assess targetability and cut probability
        if net_indel_size > 5:  # big insertion
            self.targetable = False
            self.cut_probability = 0.0
        else:
            self.compute_and_assign_cut_probability(dt)

        # update state properties
        self.cut_position = None
        self.repaired = True
        self.shift += net_indel_size


class Domain(object):
    # Each Domain may contain targets and belongs to a Genome

    def __init__(self, label, domain_start, domain_end, domain_type, genome, promoter=None):
        assert domain_type in ["orf", "promoter", "ncr"]
        self.label = label  # string
        self.domain_type = domain_type  # 'orf' or 'promoter' or 'ncr'
        self.domain_start = domain_start  # int
        self.domain_end = domain_end  # int
        if domain_type == 'orf':
            assert promoter is not None
            self.promoter = promoter  # promoter is a domain too
        self.sequence = None  # to be implemented
        self.functional = True  # bool
        self.targets = {}  # dict of Target objects and locations with labels as keys
        self.genome = genome
        genome.add_domain(self)
    
    def add_target(self, target):
        assert type(target) is Target
        assert target.domain is self
        self.targets[target.label] = target

    def update_functionality(self):
        if self.domain_type == "orf":
            if (not self.promoter.is_functional()) or (sum(target.get_shift() for target in self.targets.values()) % 3 != 0):
                self.functional = False
            else:
                self.functional = True
        elif self.domain_type == "promoter":  # TODO  how to define functional promoter
            self.functional = True
        else:  # TODO how to define functional NCR
            self.functional = True

    def genome_repair(self, label, cut_position):  # TODO FIX this method, broke it with latest PR
        location = self.target_location(label)
        location, net_indel_size, sequence = self.genome.repair_target(location, cut_position)
        self.set_location(label, location)
        # shift location of all targets to the right by net_indel_size
        for target_label in self.targets.keys():
            if self.targets[label].current_start > location:
                self.set_location(label, self.target_location(label) + net_indel_size)

        return [net_indel_size, sequence]

    def target_location(self, target_label):
        target = self.targets[target_label]
        return target.current_start

    def set_location(self, label, location):  # TODO FIX this method, broke it with latest PR
        self.targets[label][1] = location


class Genome(object):
    # Each Genome has >=1 Domains

    def __init__(self, sequence):
        self.length = len(sequence)  # int
        self.initial_genome = sequence  # string
        self.current_genome = sequence  # string
        self.repaired = True  # bool
        self.domains = {}  # dict of all domains (ORFs, promoters, NCRs)
    
    def add_domain(self, domain):
        assert type(domain) is Domain
        self.domains[domain.label] = domain

    def repair_target(self, location, cut_position, direction, sequence):  # TODO pass target instead, clean this method
        # sample from indel distribution to get left/right deletion sizes and insertion nucleotides
        # TODO: make indel() actually good
        if direction == 1:
            del_left, del_right, insert = indel()  # e.g. 0, 0, 2
        else:
            del_right, del_left, insert = indel()  # e.g. 0, 0, 2
        insert_nt = nt_rand(insert)  # fill in random sequence
        net_indel_size = insert - del_left - del_right
        left_genome = self.current_genome[0: cut_position - del_left]  # genome to left of sequence
        right_genome = self.current_genome[cut_position + del_right:]  # to right of sequence

        self.current_genome = left_genome + insert_nt + right_genome
        location = self.find_pam(location)  # fixing in case of damaged PAM
        sequence = self.current_genome[location: location+len(sequence)]
        return [location, net_indel_size, sequence]

    def find_pam(self, location):
        shift = 0
        # expands to left and right looking for nearest working PAM
        while self.current_genome[location+shift: location+shift+2] != "gg" and self.current_genome[location-shift: location-shift+2] != "gg":
            shift += 1
        if self.current_genome[location-shift: location-shift+2] == "gg":  # if nearest PAM is on left
            location -= shift  # shift location to the left
        else:  # if nearest PAM is on right
            location += shift  # shift location to the right
        return location

    def get_targets_from_genome(self):
        """Get a layered dictionary of all the targets in each domain of the genome
        Notes:
        - structure is {domain_label: dict_of_domain_targets, ...}
        """
        return {key: self.domains[key].targets for key in self.domains.keys()}

    def get_open_targets_from_genome(self):
        """Get list of unrepaired/open targets
        Notes:
        - list with format [(key_domain, key_target), ...]
        """
        open_targets = []
        target_dict = self.get_targets_from_genome()
        for key_domain in target_dict.keys():
            for key_target in target_dict[key_domain].keys():
                if not target_dict[key_domain][key_target].repaired:
                    open_targets.append((key_domain, key_target))
        return open_targets

    def get_closed_targets_from_genome(self):
        """Get list of repaired/closed targets
        Notes:
        - list with format [(key_domain, key_target), ...]
        """
        closed_targets = []
        target_dict = self.get_targets_from_genome()
        for key_domain in target_dict.keys():
            for key_target in target_dict[key_domain].keys():
                if target_dict[key_domain][key_target].repaired:
                    closed_targets.append((key_domain, key_target))
        return closed_targets

    def initialize_target_cut_probabilities(self, dt):
        """Fill in all the target cut probabilities based on dt
        Notes:
        - target cut probability initializes to None because the class doesn't naturally have access to dt
        """
        target_dict = self.get_targets_from_genome()
        for key_domain in target_dict.keys():
            for key_target in target_dict[key_domain].keys():
                target = target_dict[key_domain][key_target]
                target.compute_and_assign_cut_probability(dt)
