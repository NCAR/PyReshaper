'''
Unit tests for the Messenger class (serial only)

-----------------------
Created on May 13, 2014

@author: kpaul
'''
import unittest
from pyreshaper import messenger


class MessengerTests(unittest.TestCase):
    '''
    Messenger unit tests.
    '''

    def test_serial_init(self):
        msngr = messenger.Messenger()
        self.assertIsInstance(msngr, messenger.Messenger,
                              'Failed to create class instance')
        self.assertEqual(msngr._mpi_rank, 0,
                         'Rank is wrong after initialization')
        self.assertEqual(msngr._mpi_size, 1,
                         'Size is wrong after initialization')
        self.assertEqual(msngr._is_master, (0 == 0),
                         'Is_master is wrong after initialization')
        self.assertEqual(msngr.verbosity, 1,
                         'Verbosity is wrong after initialization')

    def test_serial_partition(self):
        msngr = messenger.Messenger()
        data = [1, 2, 3]
        p_data = msngr.partition(data)
        self.assertListEqual(data, p_data,
                         'Serial partition is wrong')

    def test_serial_is_master(self):
        msngr = messenger.Messenger()
        self.assertTrue(msngr.is_master(),
                        'Serial messenger should be master')

    def test_serial_get_rank(self):
        msngr = messenger.Messenger()
        self.assertEqual(msngr.get_rank(), 0,
                        'Serial messenger rank should be 0')

    def test_serial_get_size(self):
        msngr = messenger.Messenger()
        self.assertEqual(msngr.get_size(), 1,
                        'Serial messenger size should be 1')

    def test_serial_sum_list(self):
        msngr = messenger.Messenger()
        data = [1, 2, 3]
        msngr_sum = msngr.sum(data)
        print msngr_sum
        self.assertEqual(msngr_sum, 6,
                        'Serial messenger list sum not working')

    def test_serial_sum_dict(self):
        msngr = messenger.Messenger()
        data = {'a': 1, 'b': 2, 'c': 3}
        msngr_sum = msngr.sum(data)
        print msngr_sum
        self.assertDictEqual(msngr_sum, data,
                        'Serial messenger dict sum not working')

    def test_serial_max_list(self):
        msngr = messenger.Messenger()
        data = [1, 2, 3]
        msngr_max = msngr.sum(data)
        print msngr_max
        self.assertEqual(msngr_max, 6,
                        'Serial messenger list max not working')

    def test_serial_max_dict(self):
        msngr = messenger.Messenger()
        data = {'a': 1, 'b': 2, 'c': 3}
        msngr_max = msngr.sum(data)
        print msngr_max
        self.assertDictEqual(msngr_max, data,
                        'Serial messenger dict max not working')

    def test_serial_print_once(self):
        msngr = messenger.Messenger()
        msg = 'TEST - ONCE - SERIAL'
        msngr.print_once(msg, vlevel=0)

    def test_serial_print_all(self):
        msngr = messenger.Messenger()
        msg = 'TEST - ALL - SERIAL'
        msngr.print_all(msg, vlevel=0)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
