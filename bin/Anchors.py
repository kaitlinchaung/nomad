#!/usr/bin/env python3


class Anchor:
    def __init__(self, seq):
        self.seq = seq
        self.phase_2 = False
        self.ignorelisted = False
        self.active = False

    def is_phase_2(self):
        if self.phase_2:
            return True
        else:
            return False

    def assign_phase_2(self):
        self.phase_2 = True

    def is_ignorelisted(self):
        if self.ignorelisted:
            return True
        else:
            return False

    def ignorelist(self):
        self.ignorelisted = True

    def is_active(self):
        if self.is_active:
            return True
        else:
            return False

    def activate(self):
        self.active = True

    def deactivate(self):
        self.active = False
