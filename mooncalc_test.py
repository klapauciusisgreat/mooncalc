import mooncalc
import time
import math
import unittest


def makeUTCTime(s):
    return time.mktime(time.strptime(s, '%Y-%m-%d%Z'))-time.timezone

def toUTCString(t):
    return time.strftime('%a, %d %b %Y %H:%M:%S %Z', time.gmtime(t))

date = makeUTCTime('2013-03-05UTC')
lat = 50.5
lng = 30.5
height = 2000


class TestStringMethods(unittest.TestCase):

    def test_moon_position(self):
        moon_pos = mooncalc.get_moon_position(date, lat, lng)
        self.assertAlmostEqual(moon_pos.azimuth, -0.9783999522438226)
        self.assertAlmostEqual(moon_pos.altitude, 0.014551482243892251)
        self.assertAlmostEqual(moon_pos.distance, 364121.37256256194)

    def test_moon_illumination(self):
        moon_illum = mooncalc.get_moon_illumination(date);

        self.assertAlmostEqual(moon_illum.fraction, 0.4848068202456373)
        self.assertAlmostEqual(moon_illum.phase, 0.7548368838538762)
        self.assertAlmostEqual(moon_illum.angle, 1.6732942678578346)

    def test_get_test_moon_times(self):
        moon_times = mooncalc.get_moon_times(
            makeUTCTime('2013-03-04UTC'), lat, lng, True)
        self.assertEqual(toUTCString(moon_times.rise),
                    'Mon, 04 Mar 2013 23:54:29 UTC')
        self.assertEqual(toUTCString(moon_times.set),
                    'Mon, 04 Mar 2013 07:47:58 UTC')

if __name__ == '__main__':
    unittest.main()
