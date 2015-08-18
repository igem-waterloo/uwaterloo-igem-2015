from probabilistic import prob_cut, nt_rand, indel


# Script to store principle simulation classes and their interactions
#
# Gene/genome de-activation defined by
#   - frameshift mutation in an ORF
#   - deletion of a promoter

# TODO
# function for cut probability: foo_cut_prob(self.grna, self.target, dt, complex_concentration)
# function for cut location: foo_cut_posn()
# function for indel type (L/R del, insert) foo_indel(), foo_nt_rand(insert) etc
# if repair is big deletion, what do we fill in at the end of the target? cant be random, must come from sequence data
# need to update all genome locations with each indel or deletion
# code for deletions
# see specific TODOs / 'to be implemented' throughout the code


class Target(object):

    def __init__(self, label, grna, sequence, start, complex_concentration, direction, dt, domain):
        self.label = label  # string
        self.grna = grna  # string
        self.sequence = sequence  # string, include PAM, should be ~ 23 chars, maybe need buffer on opposite end
        self.start = start  # int
        self.total_cuts = 0  # int
        self.cut_position = None  # number of nt from PAM site, starting at 0
        self.repaired = True  # defined by open/closed
        self.targetable = True  # defined by targetable or not (PAM broken or indel size > 5)
        self.complex_concentration = complex_concentration  # conc of gRNA-cas9 complex inside nucleus
        self.shift = 0  # defined by sum of net indel sizes, used to compute frameshift if orf region
        assert direction in [1, -1]
        self.direction = direction  # 1 or -1, direction (sense) target is pointing in the genome
        self.domain = domain
        domain.add_target(self, start)
        self.cut_probability = self.compute_cut_probability(dt)

    def is_repaired(self):
        return self.repaired

    def is_targetable(self):
        return self.targetable

    def get_shift(self):
        return self.shift

    def compute_cut_probability(self, dt):  # TODO fix time dependence scope
        return prob_cut(self.grna[3:], self.sequence[3:], self.complex_concentration, dt)

    def cut(self):
        self.total_cuts += 1
        self.cut_position = 6  # 0-indexed posn of the nt right of the cut, usually 3-4nt from pam: foo_cut_posn()
        self.repaired = False

    def repair(self, dt):
        # call Genome.target_repair through Domain
        net_indel_size, self.sequence = self.domain.genome_repair(self.label, self.cut_position)

        # assess targetability and cut probability
        if net_indel_size > 5 or self.sequence[0:2] != "gg":  # big insertion or broken PAM
            self.targetable = False
            self.cut_probability = 0.0
        else:
            self.cut_probability = self.compute_cut_probability(dt)

        # update state properties
        self.cut_position = None
        self.repaired = True
        self.shift += net_indel_size


# Each Domain has >=1 targets and belongs to a Genome
class Domain(object):

    def __init__(self, label, domain_start, domain_end, domain_type, genome, promoter=None):
        assert domain_type in ["orf", "promoter", "ncr"]
        self.label = label  # string
        self.domain_type = domain_type # 'orf' or 'promoter' or 'ncr'
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
    
    def add_target(self, target, location):
        assert type(target) is Target
        assert target.domain is self
        assert type(location) is int
        self.targets[target.label] = [target, location]

    def update_functionality(self):
        if self.domain_type == "orf":
            if not self.promoter.is_functional() or sum(target.get_shift() for target in self.targets) % 3 != 0:
                self.functional = False
            else:
                self.functional = True
        elif self.domain_type == "promoter":  # TODO  how to define functional promoter
            self.functional = True
        else:  # TODO how to define functional NCR
            self.functional = True

    def genome_repair(self, label, cut_position):
        location = self.target_location(label)
        location, net_indel_size, sequence = self.genome.repair_target(location, cut_position)
        self.set_location(label, location)
        # shift location of all targets to the right by net_indel_size
        for label in self.targets:
            if self.target_location(label) > location:
                self.set_location(label, self.target_location(label)+net_indel_size)

        return [net_indel_size, sequence]

    def target_location(self, label):
        return self.targets[label][1]

    def set_location(self, label, location):
        self.targets[label][1] = location


# Each Genome has >=1 Domains
class Genome(object):

    def __init__(self, sequence):
        self.length = len(sequence)  # int
        self.initial_genome = sequence  # string
        self.current_genome = sequence  # string
        self.repaired = True  # bool
        self.domains = {}  # list of all domains (ORFs, promoters, NCRs)
    
    def add_domain(self, domain):
        assert type(domain) is Domain
        self.domains[domain.label] = domain

    def repair_target(self, location, cut_position):
        # sample from indel distribution to get left/right deletion sizes and insertion nucleotides
        # TODO: make indel() actually good
        del_left, del_right, insert = indel()  # e.g. 0, 0, 2
        insert_nt = nt_rand(insert)  # fill in random sequence
        net_indel_size = insert - del_left - del_right
        left_genome = self.current_genome[0:location+cut_position-del_left] # genome to left of sequence
        right_genome = self.current_genome[location+cut_position+del_right:] # to right of sequence

        self.current_genome = left_genome + insert_nt + right_genome
        location = self.find_pam(location)  # fixing in case of damaged PAM
        sequence = self.current_genome[location: location+23]
        return [location, net_indel_size, sequence]

    def find_pam(self, location):
        shift = 0
        # expands to left and right looking for nearest working PAM
        while self.current_genome[location+shift: location+shift+2] != "gg" and self.current_genome[location-shift: location-shift+2] != "gg":
            shift += 1
        # if nearest PAM is on left
        if self.current_genome[location-shift: location-shift+2] == "gg":
            location -= shift # shift location to the left
        else: # if nearest PAM is on right
            location += shift # shift location to the right

        return location
