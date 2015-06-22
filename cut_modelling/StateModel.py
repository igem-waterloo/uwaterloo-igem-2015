from transitions import Machine
import random

class TargetSite(object):

    states = ['cut', 'uncut']

    def __init__(self, sequence):

        self.gRNA = sequence
        self.sequence = sequence
        self.functional = True

        self.calculate_probability()

        # Initialize the state machine
        self.machine = Machine(model=self, states=TargetSite.states, initial='uncut')

        # add some transitions. We could also define these using a static list of 
        # dictionaries and then pass the list to 
        # the Machine initializer as the transitions= argument.

        self.machine.add_transition('slice', 'uncut', 'cut')

        self.machine.add_transition('repair', 'cut', 'uncut',
                                    after='indels')

    def indels(self):
        self.sequence = "abc"## self.sequence modified
        self.calculate_probability()
        self.working()

    def calculate_probability(self):
        self.cut_probability = 0.5## some function of the sequence
                               ## maybe we can also store the original
                               ## sequence (gRNA) for reference

    def working(self):
        # only checking for frameshifts
        if (len(self.sequence)-len(self.gRNA))%3 == 0:
            self.functional = True
        else:
            self.functional = False
