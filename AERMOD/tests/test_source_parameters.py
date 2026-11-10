"""Small tests of point/area allocation and AERMOD source formatting."""
from pathlib import Path
import json,tempfile,unittest,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import source_parameters as model


class SourceParametersTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[1]);self.folder=Path(self.temp.name)
        release=[dict(model_source_id='P1',source_type='POINT',weight=.25,easting_m=100.,northing_m=200.,
            elevation_m=300.,height_m=10.,temperature_K=400.,velocity_m_s=5.,diameter_m=2.5),
            dict(model_source_id='A1',source_type='AREA',weight=.75,easting_m=110.,northing_m=220.,
            elevation_m=300.,height_m=3.,x_length_m=10.,y_length_m=20.,angle_degrees=0.)]
        model.ACTIVE=(2023,'NO2','123');model.ENGINE={'WORKFLOW':'production'}
        model.PARAMETERS={'cases':{'production|2023|NO2|123':'one'},'parameters':{'one':dict(latitude=40.,longitude=-80.,elevation_m=300.,releases=release)}}
    def tearDown(self):self.temp.cleanup()
    def test_hourly_mass(self):
        file=Path(model.write_hourly(self.folder,'123',2023,{'2023-01-01':2.},{}))
        rows=[line.split() for line in file.read_text().splitlines()]
        self.assertEqual(len(rows),48)
        # Each hour: point g/s plus area g/s/m2 integrated over its 200m2 area.
        rate=float(rows[0][7])+float(rows[1][7])*200
        self.assertAlmostEqual(rate,2.)
        self.assertAlmostEqual(rate*86400,172800.)
    def test_source_geometry(self):
        def original(*args):
            file=self.folder/'model.inp'
            file.write_text('CO STARTING\nCO FINISHED\nSO STARTING\n LOCATION OLD POINT 0 0 0\n GASDEPOS OLD 0.0 0.0 0.0\nSO FINISHED\nME STARTING\nME FINISHED\n')
            return str(file)
        model.ENGINE['_write_control_file']=original
        file=Path(model.write_input('123',2023,'hourly.txt'))
        text=file.read_text()
        self.assertIn('LOCATION A1 AREA 110.00000000 220.00000000 300.00000000',text)
        self.assertIn('SRCPARAM P1 0.0 10.00000000 400.00000000 5.00000000 2.50000000',text)
        self.assertIn('SRCPARAM A1 0.0 3.00000000 10.00000000 20.00000000 0.00000000',text)
        self.assertIn('GASDEPOS A1',text)
        self.assertEqual(model.coordinates('123',{}),(40.,-80.,300.))
    def test_invalid_weights(self):
        model.PARAMETERS['parameters']['one']['releases'][0]['weight']=.1
        with self.assertRaises(ValueError):model.write_hourly(self.folder,'123',2023,{'2023-01-01':2.},{})
    def test_unknown_case(self):
        model.ACTIVE=(2024,'NO2','123')
        with self.assertRaises(ValueError):model.case_parameters()


if __name__=='__main__':unittest.main()
