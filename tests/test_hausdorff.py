from copal import hausdorff
import unittest
import numpy as np
import pandas as pd
import time

class TestHausdorff(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up for hausdorff test..')

        cls.list1 = [1,3,5,65,13,5346,7,8,7,5,4,33,6,7]
        cls.list2 = [76,5,7,6,1,423,547,78]
        cls.list3 = [1,1,1,2]
        cls.list4 = [1,1,1,-10]

        cls.array1 = np.zeros((1000,100))
        cls.array2 = np.ones((1000,100))

    def setUp(self):
        self._started_at = time.time()

    def tearDown(self):
        elapsed = time.time() - self._started_at
        print('{}: {}s'.format(self.id(),round(elapsed,2)))

    def test_hausdorff(self):
        self.assertEqual(hausdorff.hausdorff(self.list1, self.list2),
                        4799.000104188372)
        self.assertEqual(hausdorff.hausdorff(self.list3, self.list4),
                        11.045361017187261)
        
        for ix,series in enumerate(self.array1):
            self.assertEqual(hausdorff.hausdorff(series,self.array2[ix]),1)

   
    def test_new_hausdorff(self):
        self.assertEqual(hausdorff.new_hausdorff(self.list1,self.list2),
                        4799.000104188372)
        self.assertEqual(hausdorff.new_hausdorff(self.list3,self.list4),
                        11.045361017187261)

        for ix,series in enumerate(self.array1):
            self.assertEqual(hausdorff.new_hausdorff(series,self.array2[ix]),1)


if __name__ == "__main__":
    unittest.main()


