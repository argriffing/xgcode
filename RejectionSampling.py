
import unittest
import random


class RejectionSampler:

    def __init__(self, propose, accept):
        """
        @param propose: a callable that draws from the proposal distribution
        @param accept: a callable that decides if a proposal is accepted
        """
        self.propose = propose
        self.accept = accept

    def get_sample_or_none(self):
        """
        @return: a sample if accepted or None if rejected
        """
        sample = self.propose()
        if self.accept(sample):
            return sample
        else:
            return None

    def gen_samples_or_none(self, proposal_limit=None):
        """
        Each proposal yields a sample if accepted or None if rejected.
        @param proposal_limit: the number of proposals or None for no limit
        """
        proposal_count = 0
        while True:
            if proposal_limit and proposal_count >= proposal_limit:
                return
            yield self.get_sample_or_none()
            proposal_count += 1

    def gen_samples(self, acceptance_limit=None):
        """
        Each accepted proposal yields a sample.
        This function may spend a long time between yields if the acceptance rate is low.
        @param acceptance_limit: the number of accepted samples or None for no limit
        """
        acceptance_count = 0
        while True:
            if acceptance_limit and acceptance_count >= acceptance_limit:
                return
            sample = self.get_sample_or_none()
            if sample is not None:
                yield sample
                acceptance_count += 1


def propose_d6():
    """
    This is an example of a proposal function for the rejection sampler.
    It returns the result of rolling a six sided die.
    @return: an integer between 1 and 6 inclusive
    """
    return random.randint(1, 6)

def accept_d5(sample):
    """
    This is an example of an acceptance function for the rejection sampler.
    It takes an integer between 1 and 6 inclusive and returns True if it is not 6.
    @param sample: an integer between 1 and 6 inclusive
    """
    return sample < 6


class TestRejectionSampler(unittest.TestCase):

    def test(self):
        """
        Use rejection sampling to uniformly sample integers between 1 and 5 using a 6-sided die.
        """
        # create the example rejection sampler from a proposal function and an acceptance function
        sampler = RejectionSampler(propose_d6, accept_d5)
        # assert that the result of a possibly rejected sample is valid
        self.failUnless(sampler.get_sample_or_none() in (1, 2, 3, 4, 5, None))
        # after 600 proposals each valid result should have appeared at least once
        observed_set = set(sampler.gen_samples_or_none(600))
        self.failUnless(observed_set == set((1, 2, 3, 4, 5, None)))
        # after 600 accepted samples each valid result should have appeared at least once
        observed_set = set(sampler.gen_samples(600))
        self.failUnless(observed_set == set((1, 2, 3, 4, 5)))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRejectionSampler)
    unittest.TextTestRunner(verbosity=2).run(suite)

