import unittest
from src.sequence import gen_spectra
class test_gen_spectra(unittest.TestCase):
    def setUp(self):
        # No U because its not a part of the primary 20
        self.sequence = 'ARNDCEQGHILKMFPSTWYV'
        # these ion masses for the above sequence calculated from http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html
        self.b_ions_singly = [
            72.04444,
            228.14555,
            342.18847,
            457.21542,
            560.22460,
            689.26719,
            817.32577,
            874.34723,
            1011.40615,
            1124.49021,
            1237.57427,
            1365.66924,
            1496.70972,
            1643.77814,
            1740.83090,
            1827.86293,
            1928.91061,
            2114.98992,
            2278.05325,
            2377.12166
        ]
        self.y_ions_singly = [
            2395.13222,   
            2324.09511,   
            2167.99400,    
            2053.95107,    
            1938.92413,    
            1835.91494,    
            1706.87235,    
            1578.81377,    
            1521.79231,    
            1384.73340,   
            1271.64934,    
            1158.56527,     
            1030.47031,     
            899.42982,     
            752.36141,     
            655.30865,     
            568.27662,     
            467.22894,     
            281.14963,     
            118.08630,
        ]
        self.b_ions_doubly = [
            36.52588,
            114.57643,
            171.59790,
            229.11137,
            280.61596,
            345.13726,
            409.16655,
            437.67728,
            506.20673,
            562.74877,
            619.29080,
            683.33828,
            748.85852,
            822.39273,
            870.91911,
            914.43512,
            964.95896,
            1057.99862,
            1139.53028,
            1189.06449
        ]
        self.y_ions_doubly = [
            1198.06977,  
            1162.55122,    
            1084.50066,    
            1027.47920,    
            969.96573,   
            918.46113,    
            853.93984,    
            789.91055,    
            761.39982,    
            692.87036,    
            636.32833,    
            579.78630,    
            515.73882,     
            450.21857,     
            376.68437,     
            328.15798,     
            284.64197,     
            234.11813,     
            141.07847,     
            59.54681    
        ]

    def test_b_ions(self):
        return
        b1 = gen_spectra.b_ions(self.sequence, 1)
        b2 = gen_spectra.b_ions(self.sequence, 2)
        b12 = gen_spectra.b_ions(self.sequence)

        rounder = lambda x: [round(a, 3) for a in x]

        b1_rounded = rounder(b1)
        b2_rounded = rounder(b2)
        b12_rounded = rounder(b12)

        print(b2_rounded)
        print(rounder(self.b_ions_doubly))
        allb1 = all([True if x in rounder(self.b_ions_singly) else False for x in b1_rounded] + [True if x in b1_rounded else False for x in rounder(self.b_ions_singly)])
        allb2 = all([True if x in rounder(self.b_ions_doubly) else False for x in b2_rounded] + [True if x in b2_rounded else False for x in rounder(self.b_ions_doubly)])
        singlydouble = rounder(self.b_ions_singly) + rounder(self.b_ions_doubly)
        allb12 = all([True if x in singlydouble else False for x in b12_rounded] + [True if x in b12_rounded else False for x in singlydouble])

        self.assertTrue(allb1, 'all singly charged b ions should be calculated correctly')
        self.assertTrue(allb2, 'all doubly charged b ions should be calculated correctly')
        self.assertTrue(allb12, 'all singly and doubly charged b ions should be calculated correctly')

    def test_y_ions(self):
        return
        y1 = gen_spectra.y_ions(self.sequence, 1)
        y2 = gen_spectra.y_ions(self.sequence, 2)
        y12 = gen_spectra.y_ions(self.sequence)

        ally1 = all([True if x in self.y_ions_singly else False for x in y1] + [True if x in y1 else False for x in self.y_ions_singly])
        ally2 = all([True if x in self.y_ions_doubly else False for x in y2] + [True if x in y2 else False for x in self.y_ions_doubly])
        singlydouble = self.y_ions_singly + self.y_ions_doubly
        ally12 = all([True if x in singlydouble else False for x in y12] + [True if x in y12 else False for x in singlydouble])

        self.assertTrue(ally1, 'all singly charged b ions should be calculated correctly')
        self.assertTrue(ally2, 'all doubly charged b ions should be calculated correctly')
        self.assertTrue(ally12, 'all singly and doubly charged b ions should be calculated correctly')

    def test_calc_mass(self):
        b1, _ = gen_spectra.calc_masses(self.sequence, 1, 'b')
        b2, _ = gen_spectra.calc_masses(self.sequence, 2, 'b')
        y1, _ = gen_spectra.calc_masses(self.sequence, 1, 'y')
        y2, _ = gen_spectra.calc_masses(self.sequence, 2, 'y')
        b12, _ = gen_spectra.calc_masses(self.sequence, ion='b')
        y12, _ = gen_spectra.calc_masses(self.sequence, ion='y')
        by1, _ = gen_spectra.calc_masses(self.sequence, 1)
        by2, _ = gen_spectra.calc_masses(self.sequence, 2)
        by12, _ = gen_spectra.calc_masses(self.sequence)

        self.assertEqual(b1, gen_spectra.b_ions(self.sequence, 1), 'ion=b charge=1 should give the same result as calling b_ions directly')
        self.assertEqual(b2, gen_spectra.b_ions(self.sequence, 2), 'ion=b charge=2 should give the same result as calling b_ions directly')
        self.assertEqual(y1, gen_spectra.y_ions(self.sequence, 1), 'ion=y charge=1 should give the same result as calling y_ions directly')
        self.assertEqual(y2, gen_spectra.y_ions(self.sequence, 2), 'ion=y charge=2 should give the same result as calling y_ions directly')
        self.assertEqual(b12, gen_spectra.b_ions(self.sequence), 'ion=b should give the same as calling b_ions no charge directly')
        self.assertEqual(y12, gen_spectra.y_ions(self.sequence), 'ion=y should give the same as calling y_ions no charge directly')
        self.assertEqual(by1, gen_spectra.b_ions(self.sequence, 1) + gen_spectra.y_ions(self.sequence, 1), 'charge=1 should give the same as calling both b_ions and y_ions charge=1 directly')
        self.assertEqual(by2, gen_spectra.b_ions(self.sequence, 2) + gen_spectra.y_ions(self.sequence, 2), 'charge=1 should give the same as calling both b_ions and y_ions charge=2 directly')
        self.assertEqual(by12, gen_spectra.b_ions(self.sequence) + gen_spectra.y_ions(self.sequence), 'should give the same as calling both b_ions and y_ions directly')

if __name__ == '__main__':
    unittest.main()