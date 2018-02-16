from unittest import TestCase

from atmopy import compute_vars,constants,wrf_utils

class testcomputevars(TestCase):
    def test_compute_PR(self):
        pr,atts = compute_vars.compute_PR('./atmopy/tests/wrfout_test.nc')
        self.assertTrue(pr.mean()>=0.)

    def test_constants(self):
        earth_radius = constants.const.earth_radius
        self.assertTrue(earth_radius==6371000)

    def test_wrftime2date(self):
        date = wrf_utils.wrftime2date('./atmopy/tests/wrfout_test.nc')[0]
        self.assertTrue(date.year==2015)
