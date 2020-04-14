import unittest
from src.alignment import filtering
class test_filtering(unittest.TestCase):
    def setUp(self):
        self.score_key = 'score'
        self.scores = [
            {self.score_key: 0}, {self.score_key: 0}, {self.score_key: 0}, {self.score_key: 1}, {self.score_key: 2}
        ]

    def test_stddev_filter(self):
        '''
        equation is (# found in spec + streak #)/(len(reference)/2)
        '''
        sdevs2 = filtering.stddev_filter(self.scores, self.score_key, 2)
        sdevs1 = filtering.stddev_filter(self.scores, self.score_key, 1)
        sdevs0 = filtering.stddev_filter(self.scores, self.score_key, 0)

        self.assertEqual(len(sdevs2), 1, '2 standard deviations should give 1 score back')
        self.assertEqual(len(sdevs1), 2, '1 standard deviation should give 2 scores back')
        self.assertEqual(len(sdevs0), 5, '0 standard deviations should give 5 scores back')

        self.assertEqual(sdevs2, [self.scores[-1]], '2 standard deviations should give the score with value 2')
        self.assertEqual(sdevs1, self.scores[-2:], '1 standard deviation should give back scores with value 1 and 2')


if __name__ == '__main__':
    unittest.main()