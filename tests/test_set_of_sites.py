import unittest
from pyscses.set_of_sites import SetOfSites
from pyscses.site import Site
from unittest.mock import patch, Mock

class TestSetOfSites(unittest.TestCase):

    def test_set_of_sites_is_initialised(self):
        sites = tuple(Mock(spec=Site) for i in range(4))
        set_of_sites = SetOfSites(sites)
        self.assertEqual(set_of_sites.sites, sites)

    def test_set_of_sites_addition(self):
        sites1 = tuple(Mock(spec=Site) for i in range(4))
        sites2 = tuple(Mock(spec=Site) for i in range(4))
        set_of_sites_1 = SetOfSites(sites1)
        set_of_sites_2 = SetOfSites(sites2)
        combined_set_of_sites = set_of_sites_1 + set_of_sites_2
        self.assertEqual(combined_set_of_sites.sites, sites1 + sites2)

    def test_set_of_sites_sites_is_immutable(self):
        sites = tuple(Mock(spec=Site) for i in range(4))
        set_of_sites = SetOfSites(sites)
        with self.assertRaises(TypeError):
            set_of_sites.sites[2] = 'foo'

    def test_set_of_sites_stores_sites_as_a_tuple(self):
        sites = list(Mock(spec=Site) for i in range(4))
        set_of_sites = SetOfSites(sites)
        self.assertEqual(set_of_sites.sites, tuple(sites))
    
if __name__ == '__main__':
    unittest.main()
