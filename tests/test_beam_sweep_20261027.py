import io
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_beam_sweep_20261027 import import_bpp, quantile, tv_compare
from beam_sweep_array_20261027 import GzipBppStore


class BeamSweepTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir='/tmp')
        self.work=Path(self.tmp.name)
        self.db=sqlite3.connect(':memory:')
        self.db.execute('CREATE TABLE pairs(b INTEGER,i INTEGER,j INTEGER,p REAL,PRIMARY KEY(b,i,j)) WITHOUT ROWID')

    def tearDown(self):
        self.db.close();self.tmp.cleanup()

    def fixture(self,text,beam=25):
        local=self.work/str(beam);local.mkdir()
        path=local/'bpp.tsv';path.write_text(text)
        shared=self.work/(str(beam)+'-packed');shared.mkdir()
        row=GzipBppStore(shared).commit(dict(id=str(beam)),{'jobs/x/bpp.tsv':path})
        return dict(beam=beam,_result=row,_run=shared)

    def test_fixed_bases_excluded_and_unpaired_included(self):
        # Four free bases. Two change partner/unpaired mass by .5: TV .5 each.
        rec=('r','GAAAAC','(....)')
        a=self.fixture('>r\n1 6 1\n2 5 0.7\n',25)
        b=self.fixture('>r\n1 6 1\n2 5 0.2\n',400)
        qa=import_bpp(self.db,a,rec);qb=import_bpp(self.db,b,rec)
        scores=tv_compare(self.db,25,400,qa,qb,rec[2])
        self.assertAlmostEqual(scores['tv_mean'],.25)
        self.assertAlmostEqual(scores['tv_p90'],.5)
        self.assertAlmostEqual(scores['tv_max'],.5)
        self.assertEqual(scores['free_bases'],4)
        self.assertEqual(scores['max_pair_coordinates'],[2,5])

    def test_disjoint_coordinate_union(self):
        rec=('r','GCA','...')
        a=self.fixture('>r\n1 2 0.5\n',25)
        b=self.fixture('>r\n1 3 0.5\n',400)
        qa=import_bpp(self.db,a,rec);qb=import_bpp(self.db,b,rec)
        self.assertAlmostEqual(tv_compare(self.db,25,400,qa,qb,rec[2])['tv_mean'],.5)
        self.assertEqual(tv_compare(self.db,25,25,qa,qa,rec[2])['tv_max'],0)

    def test_bad_incident_mass_rejected(self):
        row=self.fixture('>r\n1 2 0.8\n1 3 0.8\n')
        with self.assertRaisesRegex(ValueError,'incident'):
            import_bpp(self.db,row,('r','GCA','...'))

    def test_missing_fixed_pair_rejected(self):
        row=self.fixture('>r\n2 5 0.5\n')
        with self.assertRaisesRegex(ValueError,'scaffold'):
            import_bpp(self.db,row,('r','GAAAAC','(....)'))

    def test_duplicate_coordinate_rejected(self):
        row=self.fixture('>r\n1 2 0.1\n1 2 0.2\n')
        with self.assertRaises(sqlite3.IntegrityError):
            import_bpp(self.db,row,('r','GCA','...'))

    def test_gzip_hash_and_raw_hash_checked(self):
        row=self.fixture('>r\n1 2 0.5\n')
        row['_result']['bpp_uncompressed_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'uncompressed BPP hash'):
            import_bpp(self.db,row,('r','GCA','...'))
        row['_result']['artifacts']['jobs/x/bpp.tsv.gz']['sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'artifact hash'):
            import_bpp(self.db,row,('r','GCA','...'))

    def test_quantile(self):
        self.assertAlmostEqual(quantile([1,0,0,0],.9),.7)
        self.assertIsNone(quantile([],.9))


if __name__=='__main__':
    unittest.main()
