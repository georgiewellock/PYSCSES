import unittest
from project.defect_at_site import Defect_at_Site
from project.site import Site
from unittest.mock import patch, Mock

class TestDefectAtSite( unittest.TestCase ):

    def setUp( self ):
        self.defect = Defect_at_Site( 'A', valence=1.0, mole_fraction=1.0, energy=1.5, site=Mock( spec=Site ) )

    def test_potential_energy( self ):
        phi = 0.5
        self.assertEqual( self.defect.potential_energy( phi ), 2.0 ) 

    def test_boltzmann_one( self ):
        phi = 0.5
        temp = 298.0
        with patch( 'project.defect_at_site.Defect_at_Site.potential_energy' ) as mock_potential_energy:
            mock_potential_energy.return_value = 2.0
            self.assertAlmostEqual( self.defect.boltzmann_one( phi=phi, temp=temp ), 1.4995936e-34, places=40 )
            mock_potential_energy.assert_called_with( phi )

if __name__ == '__main__':
    unittest.main()
