'''
Unit tests (serial only) for the TimeKeeper class

-----------------------
Created on May 31, 2014

@author: Kevin Paul <kpaul@ucar.edu>
'''
import unittest
from time import sleep
from pyreshaper import timekeeper


class TimeKeeperTests(unittest.TestCase):
    '''
    Tests for the TimeKeeper class
    '''

    def test_init(self):
        tk = timekeeper.TimeKeeper()
        self.assertIsInstance(tk, timekeeper.TimeKeeper,
            'TimeKeeper class instantiated incorrectly.')

    def test_start_stop_names(self):
        tk = timekeeper.TimeKeeper()
        name = 'Test Clock'
        wait_time = 0.05
        tk.start(name)
        sleep(wait_time)
        tk.stop(name)
        self.assertIn(name, tk._accumulated_times,
                      'Clock name not found in accumulated times dictionary')
        self.assertIn(name, tk._added_order,
                      'Clock name not found in added order list')
        self.assertIn(name, tk._start_times,
                      'Clock name not found in start times dictionary')

    def test_start_stop_values(self):
        tk = timekeeper.TimeKeeper()
        name = 'Test Clock'
        wait_time = 0.05
        ok_diff = 0.002
        tk.start(name)
        sleep(wait_time)
        tk.stop(name)
        self.assertAlmostEqual(wait_time, tk.get_time(name),
                               msg='Accumulated time seems off', delta=ok_diff)

    def test_start_stop_order_names(self):
        tk = timekeeper.TimeKeeper()
        name1 = 'Test Clock 1'
        name2 = 'Test Clock 2'
        wait_time = 0.01
        tk.start(name1)
        sleep(wait_time)
        tk.stop(name1)
        tk.start(name2)
        sleep(wait_time)
        tk.stop(name2)
        self.assertEqual(name1, tk._added_order[0],
                      'Clock name 1 not appropriately ordered')
        self.assertEqual(name2, tk._added_order[1],
                      'Clock name 2 not appropriately ordered')

    def test_start_stop_values2(self):
        tk = timekeeper.TimeKeeper()
        name1 = 'Test Clock 1'
        name2 = 'Test Clock 2'
        wait_time = 0.025
        ok_diff = 0.0025
        tk.start(name1)
        sleep(2 * wait_time)
        tk.start(name2)
        sleep(wait_time)
        tk.stop(name1)
        sleep(wait_time)
        tk.stop(name2)
        self.assertAlmostEqual(3 * wait_time, tk.get_time(name1),
                               msg='Accumulated time 1 seems off',
                               delta=ok_diff)
        self.assertAlmostEqual(2 * wait_time, tk.get_time(name2),
                               msg='Accumulated time 2 seems off',
                               delta=ok_diff)

    def test_reset_values(self):
        tk = timekeeper.TimeKeeper()
        name = 'Test Clock'
        wait_time = 0.05
        tk.start(name)
        sleep(wait_time)
        tk.stop(name)
        tk.reset(name)
        self.assertEqual(0, tk.get_time(name),
                         msg='Accumulated time seems off')

    def test_get_time(self):
        tk = timekeeper.TimeKeeper()
        name = 'Test Clock'
        wait_time = 0.05
        ok_diff = 0.0025
        tk.start(name)
        sleep(wait_time)
        tk.stop(name)
        self.assertAlmostEqual(wait_time, tk.get_time(name),
                         msg='Get time seems off',
                         delta=ok_diff)

    def test_get_all_times(self):
        tk = timekeeper.TimeKeeper()
        name1 = 'Test Clock 1'
        name2 = 'Test Clock 2'
        wait_time = 0.025
        ok_diff = 0.0025
        tk.start(name1)
        sleep(2 * wait_time)
        tk.start(name2)
        sleep(wait_time)
        tk.stop(name1)
        sleep(wait_time)
        tk.stop(name2)
        all_times = tk.get_all_times()
        expected_all_times = {name1: 3 * wait_time,
                              name2: 2 * wait_time}
        self.assertListEqual(expected_all_times.keys(), all_times.keys(),
                             'All times clock names are off')
        self.assertAlmostEqual(expected_all_times.values()[0],
                               all_times.values()[0],
                               msg='Accumulated time 1 seems off',
                               delta=ok_diff)
        self.assertAlmostEqual(expected_all_times.values()[1],
                               all_times.values()[1],
                               msg='Accumulated time 2 seems off',
                               delta=ok_diff)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
