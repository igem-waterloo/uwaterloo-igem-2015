# PAM constants
DNA_ALPHABET = "acgt"  # used for iterating over pam variants
PAM_TEMPLATE_SEQUENCE = "tggt"  # pam sequence for the progenitor file

# IO constants
PAM_TOOLS = ["Chimera", "3DNA"]
SCOREFILE_LINES = 6  # number of lines in score .txt files from dock_variants.py
CSV_HEADER = ['PAM_1', 'PAM_2', 'PAM_3', 'PAM_4', 'Tool', 'Init FA', 'Final FA', 'Init DNA', 'Final DNA',
              'Total Time', 'Dock Time']
CSV_HEADER_64 = ['PAM_1', 'PAM_2', 'PAM_3', 'Tool', 'Init FA', 'Final FA', 'Init DNA', 'Final DNA',
                 'Total Time', 'Dock Time']
