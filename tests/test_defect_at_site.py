import unittest
from project.defect_at_site import Defect_at_Site
from project.site import Site
from unittest.mock import patch, Mock

class TestDefectAtSite( unittest.TestCase ):

    def setUp( self ):
        self.defect = Defect_at_Site( 'A', valence=1.0, mole_fraction=1.0, energy=1.5, site=Mock( spec=Site ) )

    def test_standard_electrochemical_potential( self ):
        phi = 0.5
        self.assertEqual( self.defect.standard_electrochemical_potential( phi ), 2.0 ) 

    def test_boltzmann_one( self ):
        phi = 0.5
        temp = 298.0
        with patch( 'project.defect_at_site.Defect_at_Site.standard_electrochemical_potential' ) as mock_electrochemical_potential:
            mock_electrochemical_potential.return_value = 2.0
            self.assertAlmostEqual( self.defect.boltzmann_one( phi=phi, temp=temp ), 1.4995936e-34, places=40 )
            mock_electrochemical_potential.assert_called_with( phi )

if __name__ == '__main__':
    unittest.main()
