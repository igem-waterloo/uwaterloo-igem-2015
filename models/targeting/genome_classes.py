from probabilistic import prob_cut, nt_rand, indel


# Script to store principle simulation classes and their interactions
# Gene/genome de-activation defined by:
#   - frameshift mutation in an ORF
#   - deletion of a promoter
#   - deletion of a significant portion of the gene


class Target(object):

    def __init__(self, label, grna, sequence, start, complex_concentration, sense, domain):
        self.label = label  # string
        self.grna = grna  # string
        if sense == 1:
            self.sequence = sequence  # string, exclude PAM, should be ~ 20 chars
        else:
            self.sequence = self.convert_sense(sequence)
        self.original_start = start  # int, location of first target nucleotide adjacent to PAM, shouldn't change
        self.current_start = start  # int, location of first target nucleotide adjacent to PAM, changes with indels
        self.total_cuts = 0  # int, total time this target has been cut
        self.cut_position = None  # absolute genome location of cut
        self.repair_position = None # formerly cut position after repair
        self.repaired = True  # defined by open/closed
        self.targetable = True  # defined by targetable or not (PAM broken or indel size > 5)
        self.complex_concentration = complex_concentration  # conc of gRNA-cas9 complex inside nucleus
        self.shift = 0  # defined by sum of net indel sizes, used to compute frameshift if orf region
        assert sense in [1, -1]
        self.sense = sense  # 1 or -1, referring to top (explicit) or bottom (implicit) dna strand
        self.domain = domain
        domain.add_target(self)
        self.cut_probability = None

    def convert_sense(self, sequence):
        # flips characters and reverses string
        pairs = [['a', 't'], ['g', 'c']]
        converted_sequence = ""
        for nt in sequence:
            for pair in pairs:
                if nt in pair:
                    for item in pair:
                        if nt != item:
                            converted_sequence += item
        return converted_sequence[::-1]

    def is_repaired(self):
        return self.repaired

    def is_targetable(self):
        return self.targetable

    def get_shift(self):
        return self.shift

    def compute_and_assign_cut_probability(self, dt):
        grna = self.grna
        if self.sense == -1:
            grna = self.convert_sense(grna)
        self.cut_probability = prob_cut(grna, self.sequence, self.complex_concentration, dt)
        return self.cut_probability

    def cut(self):
        self.total_cuts += 1
        self.set_cut_position()
        self.repaired = False

    def set_cut_position(self):
        # posn of the nt right of the cut, usually 3-4nt from pam: foo_cut_posn()
        # hardcoded for now
        if self.sense == 1:  
            self.cut_position = self.current_start + 17
        else:
            self.cut_position = self.current_start + 3

    def repair(self, dt):
        # call Genome.target_repair through Domain
        net_indel_size = self.domain.genome.repair_target(self)
        # assess targetability and cut probability
        if abs(net_indel_size) > 5:  # big insertion
            self.targetable = False
            self.cut_probability = 0.0
        else:
            self.compute_and_assign_cut_probability(dt)

        # update state properties
        self.repair_position = self.cut_position
        self.cut_position = None
        self.repaired = True
        self.shift += net_indel_size

        # check domain functionality
        self.domain.update_functionality()


class Domain(object):
    # Each Domain may contain targets and belongs to a Genome

    def __init__(self, label, domain_start, domain_end, domain_type, genome, promoter=None):
        assert domain_type in ["orf", "promoter", "untracked"]  # note untracked isn't affected by cas9
        self.label = label  # string
        self.domain_type = domain_type  # 'orf' or 'promoter' or 'untracked'
        self.domain_start = domain_start  # int
        self.domain_end = domain_end  # int
        self.promoter = promoter
        if domain_type == 'orf':
            # assert promoter is not None
            # for now to test
            if self.promoter is not None:
                assert type(promoter) is Domain
                assert promoter.domain_type == 'promoter'
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

    def remove_target(self, target):
        assert type(target) is Target
        assert target.domain is self
        del self.targets[target.label]

    def update_functionality(self):
        if self.domain_type == "orf":
            if (not self.promoter.functional) or (sum(target.get_shift() for target in self.targets.values()) % 3 != 0):
                self.functional = False
            else:
                self.functional = True
        elif self.domain_type == "promoter":  # TODO  how to define functional promoter
            self.functional = True
        else:  # untracked domains always functional
            self.functional = True

    def target_location(self, target_label):
        target = self.targets[target_label]
        return target.current_start

    def set_location(self, target_label, location):
        self.targets[target_label].current_start = location


class Genome(object):
    # Each Genome has >=1 Domains

    def __init__(self, sequence):
        self.current_length = len(sequence)  # int
        self.initial_genome = sequence  # string
        self.current_genome = sequence  # string
        self.repaired = True  # bool
        self.domains = {}  # dict of all domains (ORFs, promoters, untracked sections)
    
    def add_domain(self, domain):
        assert type(domain) is Domain
        self.domains[domain.label] = domain

    def remove_domain(self, domain):
        assert type(domain) is Domain
        del self.domains[domain.label]

    def repair_target(self, target):
        # sample from indel distribution to get left/right deletion sizes and insertion nucleotides
        if target.sense == 1:
            del_left, del_right, insert = indel()  # e.g. 0, 0, 2
        else:
            del_right, del_left, insert = indel()  # e.g. 0, 0, 2
        insert_nt = nt_rand(insert)  # fill in random sequence
        net_indel_size = insert - del_left - del_right
        left_genome = self.current_genome[0: target.cut_position - del_left]  # genome to left of sequence
        right_genome = self.current_genome[target.cut_position + del_right:]  # to right of sequence

        new_genome = left_genome + insert_nt + right_genome
        # target.current_start = self.find_pam(target.current_start, target.sense)
        target.sequence = self.current_genome[target.current_start: target.current_start + 20]
        self.make_new_genome(len(left_genome), net_indel_size, new_genome)
        return net_indel_size

    def find_pam(self, location, sense):
        shift = 0
        # expands to left and right looking for nearest working PAM
        if sense == 1:
            while self.current_genome[location+shift+20: location+shift+23] != "gg" and self.current_genome[location-shift: location-shift+2] != "gg":
                shift += 1
            if self.current_genome[location-shift+20: location-shift+23] == "gg":  # if nearest PAM is on left
                location -= shift  # shift location to the left
            else:  # if nearest PAM is on right
                location += shift  # shift location to the right
        else:
            while self.current_genome[location+shift-3: location+shift] != "cc" and self.current_genome[location-shift: location-shift+2] != "gg":
                shift += 1
            if self.current_genome[location-shift-3: location-shift] == "cc":  # if nearest PAM is on left
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

    def make_new_genome(self, indel_location, indel_size, new_genome):
        """Re-index all domains and targets after single indel
        and set current_genome to new_genome

        This is a bit headache inducing
        """
        deleted_domains = []
        broken_targets = []
        target_dict = self.get_targets_from_genome()
        for key_domain in target_dict.keys():
            domain = self.domains[key_domain]
            # if domain starts past indel location
            if domain.domain_start > indel_location:
                # re-index it
                self.domains[key_domain].domain_start += indel_size
                # if this pulls start past indel location
                if self.domains[key_domain].domain_start < indel_location:
                    # set it to indel location
                    self.domains[key_domain].domain_start = indel_location
            # if domain ends past indel location
            if domain.domain_end > indel_location:
                # re-index it
                self.domains[key_domain].domain_end += indel_size
                # if this pulls end past indel location
                if self.domains[key_domain].domain_end < indel_location:
                    # set it to indel location
                    self.domains[key_domain].domain_end = indel_location
                # if start and end are the same (both location)
                if self.domains[key_domain].domain_start == self.domains[key_domain].domain_end:
                    # it has been deleted
                    deleted_domains.append(domain)
                else:
                    # otherwise check if it has targets to be re-indexed
                    for key_target in target_dict[key_domain].keys():
                        target = target_dict[key_domain][key_target]
                        if target.current_start > indel_location:
                            self.domains[key_domain].targets[key_target].current_start += indel_size
                            # if the target has been damaged or deleted
                            if self.domains[key_domain].targets[key_target].current_start < indel_location:
                                # deal with that later
                                broken_targets.append(self.domains[key_domain].targets[key_target])
        # remove all deleted domains
        for domain in deleted_domains:
            self.remove_domain(domain)
        # set new genome
        self.current_genome = new_genome
        self.current_length = len(new_genome)
        # self.repaired = True
        # delete or fix all broken targets
        for target in broken_targets:
            domain_label = target.domain.label
            # if whole target is deleted
            if target.current_start + 20 < indel_location:
                self.domains[domain_label].remove_target(target)
            # else it is just broken
            else:
                continue
                # find new pam, set new location, set new sequence
                # new_start = self.find_pam(target.current_start, target.sense)
                # self.domains[domain_label].targets[target.label].current_start = new_start
                # self.domains[domain_label].targets[target.label].sequence = self.current_genome[new_start:new_start+20]

    def large_deletion(self, target1, target2, dt):
        """Delete section between two open targets
        """
        assert not (target1.repaired or target2.repaired)
        target1.set_cut_position()  # make sure cut_positions of targets are up to date
        target2.set_cut_position()  # make sure cut_positions of targets are up to date
        location = min(target1.cut_position, target2.cut_position)
        middle = abs(target1.cut_position - target2.cut_position)
        if middle < self.current_length / 2:  # if middle is smaller, should delete
            new_genome = self.current_genome[0:location] + self.current_genome[location+middle:]
            self.make_new_genome(location, -middle, new_genome)
        else: # otherwise, should keep (delete beginning and end)
            # first delete beginning
            new_genome = self.current_genome[location:]
            self.make_new_genome(0, -location, new_genome)
            # then delete end
            new_genome = self.current_genome[0:middle]
            self.make_new_genome(middle, -(self.current_length - middle), new_genome)

        # keep the target on from the dleted portion if they have the same sense
        if target1.sense == target2.sense:
            if location <= target1.current_start <= location + middle:
                target_keep = target1
                target_discard = target2
            else:
                target_keep = target2
                target_discard = target1
            target_keep.sequence = target_keep.sequence[3:]
            target_keep.sequence = target_discard.sequence[0:3] + target_keep.sequence
            target_discard.domain.remove_target(target_discard)
            target_keep.cut_position = None
            target_keep.repaired = True
            target_keep.compute_and_assign_cut_probability(dt)
            # check domain functionality
            target_keep.domain.update_functionality()
        # targets have opposite sense and so they're both likely to be broken
        else:
            target1.domain.remove_target(target1)
            target2.domain.remove_target(target2)
